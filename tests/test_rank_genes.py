from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from src.phenotype.rank_genes import rank_genes


class RankGenesTest(unittest.TestCase):
    def test_gene_matching_multiple_groups_ranks_first(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ontology = root / "hp.obo"
            annotations = root / "genes.tsv"
            config = root / "config.yaml"

            ontology.write_text(
                """
[Term]
id: HP:0000001
name: Root

[Term]
id: HP:0000002
name: Cancer
is_a: HP:0000001 ! Root

[Term]
id: HP:0000003
name: Growth
is_a: HP:0000001 ! Root
""".strip(),
                encoding="utf-8",
            )
            annotations.write_text(
                "1\tMULTI\tHP:0000002\tCancer\n"
                "1\tMULTI\tHP:0000003\tGrowth\n"
                "2\tSINGLE\tHP:0000002\tCancer\n",
                encoding="utf-8",
            )
            config.write_text(
                """
hpo:
  - id: HP:0000002
    label: Cancer
    group: cancer
    weight: 1
  - id: HP:0000003
    label: Growth
    group: growth
    weight: 1
""".strip(),
                encoding="utf-8",
            )

            ranked = rank_genes(ontology, annotations, config)

        self.assertEqual(ranked[0]["gene"], "MULTI")
        self.assertEqual(ranked[0]["group_coverage"], 1.0)


if __name__ == "__main__":
    unittest.main()
