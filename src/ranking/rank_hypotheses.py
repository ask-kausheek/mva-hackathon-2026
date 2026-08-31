from __future__ import annotations

import argparse
import csv
from pathlib import Path

INHERITANCE_TERMS = {
    "HP:0000006": "AD",
    "HP:0000007": "AR",
    "HP:0001417": "X-linked",
    "HP:0001427": "mitochondrial",
}


def load_inheritance(path: Path) -> dict[str, set[str]]:
    inheritance: dict[str, set[str]] = {}
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) < 3:
                continue
            gene = columns[1]
            for value in columns:
                mode = INHERITANCE_TERMS.get(value)
                if mode:
                    inheritance.setdefault(gene, set()).add(mode)
    return inheritance


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_variant_key(key: str) -> tuple[str, str, str, str]:
    chrom, pos, ref, alt = key.split(":", 3)
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    return chrom, pos, ref, alt


def write_submission(
    path: Path,
    hypotheses: list[dict[str, object]],
    proband_id: str,
) -> None:
    columns = [
        "proband_id",
        "chrom_1",
        "pos_1",
        "ref_1",
        "alt_1",
        "chrom_2",
        "pos_2",
        "ref_2",
        "alt_2",
        "epcr",
        "finding_type",
        "notes",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for hypothesis in hypotheses:
            chrom_1, pos_1, ref_1, alt_1 = parse_variant_key(
                str(hypothesis["variant_1"])
            )
            variant_2 = str(hypothesis["variant_2"])
            if variant_2:
                chrom_2, pos_2, ref_2, alt_2 = parse_variant_key(variant_2)
            else:
                chrom_2 = pos_2 = ref_2 = alt_2 = ""

            writer.writerow(
                {
                    "proband_id": proband_id,
                    "chrom_1": chrom_1,
                    "pos_1": pos_1,
                    "ref_1": ref_1,
                    "alt_1": alt_1,
                    "chrom_2": chrom_2,
                    "pos_2": pos_2,
                    "ref_2": ref_2,
                    "alt_2": alt_2,
                    "epcr": hypothesis["epcr"],
                    "finding_type": hypothesis["finding_type"],
                    "notes": (
                        f"{hypothesis['gene']} {hypothesis['hypothesis_type']}; "
                        f"{str(hypothesis['rationale']).replace(',', '; ')}"
                    ),
                }
            )


def single_hypotheses(
    rows: list[dict[str, str]],
    inheritance: dict[str, set[str]],
) -> list[dict[str, object]]:
    hypotheses: list[dict[str, object]] = []

    for row in rows:
        gene = row["gene"]
        if gene.startswith("HLA-"):
            continue
        clin_sig = row.get("clin_sig", "").lower()
        if "benign" in clin_sig:
            continue
        modes = inheritance.get(gene, set())
        genotype = row["gt"]
        base_score = float(row["single_score"])
        phenotype_rank = int(row["phenotype_rank"])
        adjustment = 0.0
        reasons: list[str] = []
        finding_type = "primary"
        is_pathogenic = (
            ("pathogenic" in clin_sig or "likely_pathogenic" in clin_sig)
            and "conflict" not in clin_sig
        )
        is_vus = "uncertain_significance" in clin_sig or "conflict" in clin_sig

        if genotype in {"1/1", "1|1"}:
            hypothesis_type = "homozygous"
            if "AR" in modes:
                adjustment += 0.12
                reasons.append("homozygous_AR_compatible")
            else:
                adjustment -= 0.12
                reasons.append("homozygous_without_AR_support")
            if phenotype_rank > 500:
                adjustment -= 0.20
                finding_type = "secondary"
                reasons.append("weak_case_phenotype_rank")
        else:
            hypothesis_type = "single_heterozygous"
            if modes == {"AR"}:
                adjustment -= 0.25
                reasons.append("single_allele_in_AR_gene")
            elif is_pathogenic:
                if phenotype_rank <= 250:
                    adjustment += 0.05
                    reasons.append("pathogenic_with_strong_phenotype_rank")
                elif phenotype_rank <= 500:
                    adjustment -= 0.05
                    reasons.append("pathogenic_with_moderate_phenotype_rank")
                else:
                    adjustment -= 0.25
                    finding_type = "secondary"
                    reasons.append("pathogenic_but_weak_case_phenotype")
            elif is_vus:
                if phenotype_rank <= 100:
                    adjustment -= 0.05
                    reasons.append("VUS_with_strong_phenotype_rank")
                elif phenotype_rank <= 250:
                    adjustment -= 0.15
                    reasons.append("VUS_with_moderate_phenotype_rank")
                else:
                    adjustment -= 0.30
                    finding_type = "secondary"
                    reasons.append("VUS_with_weak_case_phenotype")
            else:
                if phenotype_rank <= 100:
                    adjustment -= 0.05
                    reasons.append("novel_with_strong_phenotype_rank")
                elif phenotype_rank <= 250:
                    adjustment -= 0.20
                    reasons.append("novel_with_moderate_phenotype_rank")
                else:
                    adjustment -= 0.35
                    finding_type = "secondary"
                    reasons.append("novel_with_weak_case_phenotype")

        hypotheses.append(
            {
                "hypothesis_type": hypothesis_type,
                "gene": gene,
                "variant_1": row["variant_key"],
                "variant_2": "",
                "genotype": genotype,
                "inheritance": ",".join(sorted(modes)) or "unknown",
                "phenotype_rank": phenotype_rank,
                "base_score": base_score,
                "adjustment": adjustment,
                "hypothesis_score": max(0.0, min(1.0, base_score + adjustment)),
                "finding_type": finding_type,
                "rationale": ",".join(reasons),
            }
        )

    return hypotheses


def pair_hypotheses(
    rows: list[dict[str, str]],
    inheritance: dict[str, set[str]],
) -> list[dict[str, object]]:
    hypotheses: list[dict[str, object]] = []

    for row in rows:
        gene = row["gene"]
        modes = inheritance.get(gene, set())
        base_score = float(row["pair_score"])
        phenotype_rank = int(row["phenotype_rank"])
        adjustment = 0.15
        reasons = ["two_allele_model"]

        if "AR" in modes:
            adjustment += 0.10
            reasons.append("pair_AR_compatible")
        else:
            adjustment -= 0.15
            reasons.append("pair_without_AR_support")

        if row["phase_status"] == "trans_supported":
            adjustment += 0.05
            reasons.append("trans_supported")
        elif row["phase_status"] == "unknown":
            reasons.append("phase_unknown")

        finding_type = "primary"
        if phenotype_rank > 500:
            adjustment -= 0.20
            finding_type = "secondary"
            reasons.append("weak_case_phenotype_rank")

        hypotheses.append(
            {
                "hypothesis_type": "compound_heterozygous",
                "gene": gene,
                "variant_1": row["variant_1"],
                "variant_2": row["variant_2"],
                "genotype": "0/1 + 0/1",
                "inheritance": ",".join(sorted(modes)) or "unknown",
                "phenotype_rank": phenotype_rank,
                "base_score": base_score,
                "adjustment": adjustment,
                "hypothesis_score": max(0.0, min(1.0, base_score + adjustment)),
                "finding_type": finding_type,
                "rationale": ",".join(reasons),
            }
        )

    return hypotheses


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--singles", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--gene-annotations", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--submission-output",
        type=Path,
        default=Path("submission/ask-kausheek_hpo-vep-rankfusion.csv"),
    )
    parser.add_argument("--proband-id", default="PROBAND01")
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()

    inheritance = load_inheritance(args.gene_annotations)
    pair_rows = load_rows(args.pairs)
    pair_variants = {
        row[key]
        for row in pair_rows
        for key in ("variant_1", "variant_2")
    }
    hypotheses = [
        row
        for row in single_hypotheses(load_rows(args.singles), inheritance)
        if row["variant_1"] not in pair_variants
    ]
    hypotheses.extend(pair_hypotheses(pair_rows, inheritance))
    hypotheses.sort(
        key=lambda row: (
            -float(row["hypothesis_score"]),
            0 if row["hypothesis_type"] == "compound_heterozygous" else 1,
            int(row["phenotype_rank"]),
        )
    )

    selected: list[dict[str, object]] = []
    seen_single_genes: set[str] = set()
    for hypothesis in hypotheses:
        if hypothesis["hypothesis_type"] == "single_heterozygous":
            gene = str(hypothesis["gene"])
            if gene in seen_single_genes:
                continue
            seen_single_genes.add(gene)
        selected.append(hypothesis)
        if len(selected) == args.top:
            break

    scores = [float(row["hypothesis_score"]) for row in selected]
    minimum_score = min(scores)
    maximum_score = max(scores)
    for index, row in enumerate(selected):
        if maximum_score == minimum_score:
            epcr = 1.0 - (0.95 * index / max(len(selected) - 1, 1))
        else:
            epcr = 0.05 + 0.95 * (
                (float(row["hypothesis_score"]) - minimum_score)
                / (maximum_score - minimum_score)
            )
        row["epcr"] = round(epcr, 3)

    for index in range(1, len(selected)):
        if float(selected[index]["epcr"]) >= float(selected[index - 1]["epcr"]):
            selected[index]["epcr"] = round(
                max(0.001, float(selected[index - 1]["epcr"]) - 0.001),
                3,
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    columns = ["rank", *selected[0].keys()]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for rank, hypothesis in enumerate(selected, start=1):
            writer.writerow({"rank": rank, **hypothesis})

    write_submission(args.submission_output, selected, args.proband_id)

    for rank, row in enumerate(selected, start=1):
        print(
            rank,
            row["hypothesis_type"],
            row["gene"],
            row["variant_1"],
            row["variant_2"] or "-",
            f"score={row['hypothesis_score']:.4f}",
            row["finding_type"],
            row["rationale"],
            sep="\t",
        )
    print(f"\nSubmission CSV: {args.submission_output}")


if __name__ == "__main__":
    main()
