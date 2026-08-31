from __future__ import annotations

import argparse
import math
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import yaml

HPO_PATTERN = re.compile(r"^HP:\d{7}$")


def parse_obo(path: Path) -> tuple[dict[str, str], dict[str, set[str]], dict[str, str]]:
    names: dict[str, str] = {}
    parents: dict[str, set[str]] = defaultdict(set)
    aliases: dict[str, str] = {}
    current: dict[str, object] | None = None

    def commit(term: dict[str, object] | None) -> None:
        if not term or term.get("obsolete") or "id" not in term:
            return
        term_id = str(term["id"])
        names[term_id] = str(term.get("name", term_id))
        parents[term_id].update(term.get("parents", set()))
        for alias in term.get("aliases", set()):
            aliases[str(alias)] = term_id

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            if line == "[Term]":
                commit(current)
                current = {"parents": set(), "aliases": set(), "obsolete": False}
            elif line.startswith("["):
                commit(current)
                current = None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("alt_id: "):
                    current["aliases"].add(line[8:])
                elif line.startswith("is_a: "):
                    current["parents"].add(line[6:].split(" ! ", 1)[0])
                elif line == "is_obsolete: true":
                    current["obsolete"] = True
        commit(current)

    return names, dict(parents), aliases


def parse_gene_annotations(path: Path, aliases: dict[str, str]) -> dict[str, set[str]]:
    annotations: dict[str, set[str]] = defaultdict(set)

    with path.open(encoding="utf-8-sig") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")
            hpo_indexes = [i for i, value in enumerate(columns) if HPO_PATTERN.match(value)]
            if not hpo_indexes:
                continue

            hpo_index = hpo_indexes[0]
            hpo_id = aliases.get(columns[hpo_index], columns[hpo_index])

            if hpo_index >= 2:
                symbol = columns[1].strip()
            else:
                continue

            if symbol and symbol.lower() not in {"gene_symbol", "entrez-gene-symbol"}:
                annotations[symbol].add(hpo_id)

    if not annotations:
        raise ValueError(f"No gene-HPO annotations parsed from {path}")
    return dict(annotations)


def rank_genes(
    ontology_path: Path,
    annotations_path: Path,
    config_path: Path,
) -> list[dict[str, object]]:
    names, parents, aliases = parse_obo(ontology_path)
    gene_terms = parse_gene_annotations(annotations_path, aliases)

    with config_path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)

    query = [
        {
            **item,
            "id": aliases.get(item["id"], item["id"]),
            "weight": float(item.get("weight", 1.0)),
            "group": item.get("group", "other"),
        }
        for item in config["hpo"]
    ]

    @lru_cache(maxsize=None)
    def ancestors(term_id: str) -> frozenset[str]:
        result = {term_id}
        for parent in parents.get(term_id, set()):
            result.update(ancestors(parent))
        return frozenset(result)

    propagated_by_gene: dict[str, set[str]] = {}
    term_gene_counts: dict[str, int] = defaultdict(int)
    for symbol, direct_terms in gene_terms.items():
        propagated = set().union(*(ancestors(term) for term in direct_terms))
        propagated_by_gene[symbol] = propagated
        for term in propagated:
            term_gene_counts[term] += 1

    gene_count = len(gene_terms)

    def information_content(term_id: str) -> float:
        return -math.log((term_gene_counts.get(term_id, 0) + 1) / (gene_count + 1))

    def lin_similarity(left: str, right: str) -> float:
        common = ancestors(left) & ancestors(right)
        if not common:
            return 0.0
        mica = max(information_content(term) for term in common)
        denominator = information_content(left) + information_content(right)
        return 2 * mica / denominator if denominator else 0.0

    results: list[dict[str, object]] = []
    all_groups = {item["group"] for item in query}

    for symbol, direct_terms in gene_terms.items():
        weighted_similarity = 0.0
        total_weight = 0.0
        covered_groups: set[str] = set()
        matches: list[str] = []

        for feature in query:
            best_term, best_similarity = max(
                (
                    (gene_term, lin_similarity(feature["id"], gene_term))
                    for gene_term in direct_terms
                ),
                key=lambda item: item[1],
            )
            weight = feature["weight"]
            weighted_similarity += weight * best_similarity
            total_weight += weight

            if best_similarity >= 0.5:
                covered_groups.add(feature["group"])
            if best_similarity >= 0.25:
                matches.append(
                    f"{feature['id']}->{best_term}:{best_similarity:.3f}"
                )

        semantic_score = weighted_similarity / total_weight
        coverage = len(covered_groups) / len(all_groups)
        final_score = 0.85 * semantic_score + 0.15 * coverage
        results.append(
            {
                "gene": symbol,
                "phenotype_score": final_score,
                "semantic_score": semantic_score,
                "group_coverage": coverage,
                "matched_groups": ",".join(sorted(covered_groups)),
                "matches": ";".join(matches),
                "annotation_count": len(direct_terms),
            }
        )

    return sorted(
        results,
        key=lambda item: (
            -float(item["phenotype_score"]),
            -float(item["group_coverage"]),
            str(item["gene"]),
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ontology", required=True, type=Path)
    parser.add_argument("--annotations", required=True, type=Path)
    parser.add_argument("--config", type=Path, default=Path("configs/proband.yaml"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--top", type=int, default=500)
    args = parser.parse_args()

    ranked = rank_genes(args.ontology, args.annotations, args.config)[: args.top]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "rank",
        "gene",
        "phenotype_score",
        "semantic_score",
        "group_coverage",
        "matched_groups",
        "annotation_count",
        "matches",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for rank, item in enumerate(ranked, start=1):
            row = {
                **item,
                "rank": rank,
                "phenotype_score": f"{item['phenotype_score']:.6f}",
                "semantic_score": f"{item['semantic_score']:.6f}",
                "group_coverage": f"{item['group_coverage']:.6f}",
            }
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")


if __name__ == "__main__":
    main()
