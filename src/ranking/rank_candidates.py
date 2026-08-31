from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import re
from pathlib import Path

IMPACT_SCORE = {
    "HIGH": 1.0,
    "MODERATE": 0.7,
    "LOW": 0.25,
    "MODIFIER": 0.05,
}

CONSEQUENCE_SCORE = {
    "transcript_ablation": 1.0,
    "splice_acceptor_variant": 1.0,
    "splice_donor_variant": 1.0,
    "stop_gained": 1.0,
    "frameshift_variant": 1.0,
    "stop_lost": 0.9,
    "start_lost": 0.9,
    "transcript_amplification": 0.8,
    "inframe_insertion": 0.55,
    "inframe_deletion": 0.55,
    "missense_variant": 0.7,
    "protein_altering_variant": 0.6,
    "splice_region_variant": 0.4,
    "synonymous_variant": 0.1,
}


def optional_float(value: str | None) -> float | None:
    if value in (None, "", "."):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def optional_int(value: str | None) -> int | None:
    if value in (None, "", "."):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def functional_score(row: dict[str, str]) -> float:
    consequences = set((row.get("consequence") or "").split("&"))
    consequence_score = max(
        (CONSEQUENCE_SCORE.get(item, 0.05) for item in consequences),
        default=0.05,
    )
    return max(IMPACT_SCORE.get(row.get("impact", ""), 0.0), consequence_score)


def rarity_score(max_af: float | None) -> float:
    if max_af is None:
        return 0.8
    if max_af <= 0.00001:
        return 1.0
    if max_af <= 0.0001:
        return 0.9
    if max_af <= 0.001:
        return 0.75
    if max_af <= 0.01:
        return 0.45
    if max_af <= 0.05:
        return 0.15
    return 0.0


def clinical_score(value: str | None) -> tuple[float, bool]:
    tokens = {
        token.strip().lower()
        for token in re.split(r"[,&]", value or "")
        if token.strip()
    }
    has_pathogenic = bool({"pathogenic", "likely_pathogenic"} & tokens)
    has_benign = bool({"benign", "likely_benign", "benign/likely_benign"} & tokens)
    has_conflict = any("conflict" in token for token in tokens)

    if has_conflict or (has_pathogenic and has_benign):
        return 0.2, False
    if "pathogenic" in tokens:
        return 1.0, True
    if "likely_pathogenic" in tokens:
        return 0.9, True
    if "benign" in tokens or "likely_benign" in tokens:
        return -0.5, False
    if "uncertain_significance" in tokens:
        return 0.25, False
    return 0.0, False


def technical_features(row: dict[str, str]) -> tuple[float, int, float | None]:
    dp = optional_int(row.get("dp")) or 0
    gq = optional_int(row.get("gq")) or 0
    ad_parts = [
        optional_int(value) or 0
        for value in (row.get("ad") or "").split(",")
        if value not in ("", ".")
    ]
    vaf = ad_parts[1] / sum(ad_parts) if len(ad_parts) >= 2 and sum(ad_parts) else None
    score = (
        0.35 * min(gq / 99, 1)
        + 0.30 * min(dp / 30, 1)
        + 0.20 * (1 if row.get("filter") == "PASS" else 0)
        + 0.15 * (1 if vaf is not None and 0.20 <= vaf <= 0.80 else 0)
    )
    return score, dp, vaf


def load_phenotype_scores(path: Path) -> dict[str, tuple[float, int]]:
    result: dict[str, tuple[float, int]] = {}
    with path.open(encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["gene"]] = (float(row["phenotype_score"]), int(row["rank"]))
    return result


def annotation_priority(row: dict[str, str]) -> tuple[float, int, int]:
    return (
        functional_score(row),
        1 if row.get("mane") not in (None, "", ".") else 0,
        1 if row.get("canonical") == "YES" else 0,
    )


def load_candidates(
    annotation_path: Path,
    phenotype_scores: dict[str, tuple[float, int]],
) -> list[dict[str, object]]:
    best_annotations: dict[tuple[str, str, str, str, str], dict[str, str]] = {}

    with gzip.open(annotation_path, "rt", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row.get("gene")
            if gene in (None, "", ".") or gene not in phenotype_scores:
                continue
            key = (row["chrom"], row["pos"], row["ref"], row["alt"], gene)
            current = best_annotations.get(key)
            if current is None or annotation_priority(row) > annotation_priority(current):
                best_annotations[key] = row

    candidates: list[dict[str, object]] = []
    for row in best_annotations.values():
        max_af = optional_float(row.get("max_af"))
        function = functional_score(row)
        clinical, pathogenic = clinical_score(row.get("clin_sig"))

        # Broad candidate generation: functional/splice or known pathogenic,
        # with a permissive 1% population threshold. Missing AF is retained.
        if clinical < 0:
            continue
        if function < 0.4 and not pathogenic:
            continue
        if max_af is not None and max_af > 0.01:
            continue

        phenotype, phenotype_rank = phenotype_scores[row["gene"]]
        rarity = rarity_score(max_af)
        technical, dp, vaf = technical_features(row)
        single_score = (
            0.30 * phenotype
            + 0.25 * function
            + 0.20 * rarity
            + 0.15 * clinical
            + 0.10 * technical
        )

        candidates.append(
            {
                **row,
                "variant_key": f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}",
                "phenotype_rank": phenotype_rank,
                "phenotype_score": phenotype,
                "functional_score": function,
                "rarity_score": rarity,
                "clinical_score": clinical,
                "technical_score": technical,
                "max_af_numeric": max_af,
                "depth": dp,
                "vaf": vaf,
                "single_score": single_score,
            }
        )

    return sorted(candidates, key=lambda row: -float(row["single_score"]))


def phase_status(left: dict[str, object], right: dict[str, object]) -> str:
    left_pid = left.get("pid")
    right_pid = right.get("pid")
    left_pgt = left.get("pgt")
    right_pgt = right.get("pgt")
    if left_pid in (None, "", ".") or left_pid != right_pid:
        return "unknown"
    if left_pgt in (None, "", ".") or right_pgt in (None, "", "."):
        return "unknown"
    return "trans_supported" if left_pgt != right_pgt else "cis_supported"


def build_pairs(candidates: list[dict[str, object]]) -> list[dict[str, object]]:
    by_gene: dict[str, list[dict[str, object]]] = {}
    for candidate in candidates:
        if candidate["gt"] not in {"0/1", "1/0", "0|1", "1|0"}:
            continue
        by_gene.setdefault(str(candidate["gene"]), []).append(candidate)

    pairs: list[dict[str, object]] = []
    for gene, gene_candidates in by_gene.items():
        unique: dict[str, dict[str, object]] = {}
        for candidate in gene_candidates:
            unique[str(candidate["variant_key"])] = candidate
        shortlist = sorted(
            unique.values(),
            key=lambda row: -float(row["single_score"]),
        )[:25]

        for left, right in itertools.combinations(shortlist, 2):
            phase = phase_status(left, right)
            if phase == "cis_supported":
                continue

            distance = abs(int(left["pos"]) - int(right["pos"]))
            if phase == "unknown" and distance < 500:
                # Nearby unphased calls frequently represent one complex
                # haplotype or normalization fragments, not independent alleles.
                continue

            # Pair hypotheses require stronger rarity than single-variant
            # retrieval. Known pathogenic alleles remain eligible regardless
            # of population-frequency availability.
            left_pair_eligible = (
                left["max_af_numeric"] is None
                or float(left["max_af_numeric"]) <= 0.001
                or float(left["clinical_score"]) >= 0.9
            )
            right_pair_eligible = (
                right["max_af_numeric"] is None
                or float(right["max_af_numeric"]) <= 0.001
                or float(right["clinical_score"]) >= 0.9
            )
            if not (left_pair_eligible and right_pair_eligible):
                continue

            if phase == "trans_supported":
                phase_adjustment = 0.10
            else:
                phase_adjustment = 0.0

            pair_score = (
                0.35 * (float(left["single_score"]) + float(right["single_score"])) / 2
                + 0.20
                * min(float(left["functional_score"]), float(right["functional_score"]))
                + 0.20 * min(float(left["rarity_score"]), float(right["rarity_score"]))
                + 0.15 * float(left["phenotype_score"])
                + 0.10
                * min(float(left["technical_score"]), float(right["technical_score"]))
                + phase_adjustment
            )
            pairs.append(
                {
                    "gene": gene,
                    "variant_1": left["variant_key"],
                    "variant_2": right["variant_key"],
                    "consequence_1": left["consequence"],
                    "consequence_2": right["consequence"],
                    "hgvsp_1": left["hgvsp"],
                    "hgvsp_2": right["hgvsp"],
                    "max_af_1": left["max_af_numeric"],
                    "max_af_2": right["max_af_numeric"],
                    "clin_sig_1": left["clin_sig"],
                    "clin_sig_2": right["clin_sig"],
                    "phase_status": phase,
                    "phenotype_rank": left["phenotype_rank"],
                    "pair_score": pair_score,
                }
            )

    return sorted(pairs, key=lambda row: -float(row["pair_score"]))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--phenotype-ranks", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    phenotype_scores = load_phenotype_scores(args.phenotype_ranks)
    candidates = load_candidates(args.annotations, phenotype_scores)
    pairs = build_pairs(candidates)

    write_tsv(args.output_dir / "ranked_single_candidates.tsv", candidates)
    write_tsv(args.output_dir / "ranked_pairs.tsv", pairs)

    print(f"Candidate variant-gene rows: {len(candidates)}")
    print(f"Candidate heterozygous pairs: {len(pairs)}")
    print("\nTop single candidates")
    for rank, row in enumerate(candidates[:20], start=1):
        print(
            rank,
            row["variant_key"],
            row["gene"],
            row["consequence"],
            f"score={row['single_score']:.4f}",
            f"max_af={row['max_af_numeric']}",
            row["clin_sig"],
            sep="\t",
        )
    print("\nTop pairs")
    for rank, row in enumerate(pairs[:20], start=1):
        print(
            rank,
            row["gene"],
            row["variant_1"],
            row["variant_2"],
            f"score={row['pair_score']:.4f}",
            row["phase_status"],
            sep="\t",
        )


if __name__ == "__main__":
    main()
