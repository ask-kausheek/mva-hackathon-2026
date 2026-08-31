#!/usr/bin/env bash
set -euo pipefail

python -m src.phenotype.rank_genes \
  --ontology resources/hpo/hp.obo \
  --annotations resources/hpo/genes_to_phenotype.txt \
  --config configs/proband.yaml \
  --output data/phenotype/ranked_genes.tsv \
  --top "${1:-500}"

cut -f1-7 data/phenotype/ranked_genes.tsv | head -21
