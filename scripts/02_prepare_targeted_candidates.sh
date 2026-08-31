#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 INPUT.vcf.gz REFERENCE.fa GENES.txt GENES.gtf.gz OUTPUT_DIR" >&2
  exit 2
fi

input=$1
reference=$2
genes_file=$3
gtf=$4
output_dir=$5

mkdir -p "$output_dir"
regions="$output_dir/gene_regions.bed"
subset="$output_dir/targeted_candidates.vcf.gz"

gzip -cd "$gtf" |
  awk -v genes="$genes_file" '
    BEGIN {
      FS = "\t"
      OFS = "\t"
      while ((getline line < genes) > 0) {
        if (line !~ /^#/ && line != "") wanted[line] = 1
      }
      close(genes)
    }
    $3 == "gene" {
      name = ""
      count = split($9, attributes, ";")
      for (i = 1; i <= count; i++) {
        if (attributes[i] ~ /gene_name "/) {
          name = attributes[i]
          sub(/^[[:space:]]*gene_name "/, "", name)
          sub(/"$/, "", name)
        }
      }
      if (name in wanted) print $1, $4 - 1, $5, name
    }
  ' |
  sort -k1,1V -k2,2n >"$regions"

if [[ ! -s "$regions" ]]; then
  echo "No requested genes were found in the GTF; $regions is empty." >&2
  exit 1
fi

bcftools view -R "$regions" "$input" -Ou |
  bcftools norm -f "$reference" -m -any --keep-sum AD -Ou |
  bcftools view \
    -i 'GT="het" && FMT/DP>=8 && FMT/GQ>=20 && FMT/AD[0:1]>=3' \
    -Oz -o "$subset"

tabix -f -p vcf "$subset"
bcftools stats "$subset" >"$output_dir/targeted_candidates.stats.txt"

echo "Created $regions"
echo "Created $subset"
