#!/usr/bin/env bash
set -euo pipefail

annotations=${1:-data/annotation/hpo_top2000/hpo_candidates.vep.tsv.gz}
output_dir=${2:-data/ranking/top2000}

python -m src.ranking.rank_candidates \
  --annotations "$annotations" \
  --phenotype-ranks data/phenotype/ranked_genes.tsv \
  --output-dir "$output_dir"
