#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 || $# -gt 7 ]]; then
  cat >&2 <<'EOF'
Usage:
  05_prepare_hpo_candidates.sh RANKED_GENES.tsv GENES.gtf.gz INPUT.vcf.gz \
    REFERENCE.fa OUTPUT_DIR [TOP_GENES] [SPLICE_FLANK]
EOF
  exit 2
fi

ranked_genes=$1
gtf=$2
input=$3
reference=$4
output_dir=$5
top_genes=${6:-500}
splice_flank=${7:-20}

for command in bcftools bgzip tabix gzip awk sort head tail cut sed grep; do
  command -v "$command" >/dev/null || {
    echo "Missing required command: $command" >&2
    exit 1
  }
done

for file in "$ranked_genes" "$gtf" "$input" "$reference"; do
  [[ -f "$file" ]] || {
    echo "Required file not found: $file" >&2
    exit 1
  }
done

mkdir -p "$output_dir"
gene_list="$output_dir/hpo_top${top_genes}_genes.txt"
raw_regions="$output_dir/hpo_top${top_genes}_exons.raw.bed"
regions="$output_dir/hpo_top${top_genes}_coding_splice.bed"
candidates="$output_dir/hpo_top${top_genes}_candidates.vcf.gz"

tail -n +2 "$ranked_genes" |
  head -n "$top_genes" |
  cut -f2 |
  sed 's/\r$//' >"$gene_list"

gzip -cd "$gtf" |
  awk -v genes="$gene_list" -v flank="$splice_flank" '
    BEGIN {
      FS = "\t"
      OFS = "\t"
      while ((getline line < genes) > 0) {
        if (line != "") wanted[line] = 1
      }
      close(genes)
    }
    $3 == "exon" {
      name = ""
      count = split($9, attributes, ";")
      for (i = 1; i <= count; i++) {
        if (attributes[i] ~ /gene_name "/) {
          name = attributes[i]
          sub(/^[[:space:]]*gene_name "/, "", name)
          sub(/"$/, "", name)
        }
      }
      if (name in wanted) {
        start = $4 - 1 - flank
        if (start < 0) start = 0
        print $1, start, $5 + flank
      }
    }
  ' |
  sort -k1,1V -k2,2n -k3,3n >"$raw_regions"

if [[ ! -s "$raw_regions" ]]; then
  echo "No exon regions were found for the ranked genes." >&2
  exit 1
fi

# Merge overlapping transcript exons without requiring bedtools.
awk '
  BEGIN { OFS = "\t" }
  NR == 1 {
    chrom = $1
    start = $2
    end = $3
    next
  }
  $1 == chrom && $2 <= end {
    if ($3 > end) end = $3
    next
  }
  {
    print chrom, start, end
    chrom = $1
    start = $2
    end = $3
  }
  END {
    if (NR > 0) print chrom, start, end
  }
' "$raw_regions" >"$regions"

# Keep all technically credible non-reference genotypes. Genotype class is
# retained for later heterozygous, homozygous and low-VAF hypothesis scoring.
bcftools view -R "$regions" "$input" -Ou |
  bcftools norm -f "$reference" -m -any --keep-sum AD -Ou |
  bcftools view \
    -i 'GT~"[1-9]" && FMT/DP>=8 && FMT/GQ>=20 && FMT/AD[0:1]>=3' \
    -Oz -o "$candidates"

tabix -f -p vcf "$candidates"
bcftools stats "$candidates" >"$output_dir/hpo_top${top_genes}_candidates.stats.txt"

echo "Genes:      $gene_list"
echo "Regions:    $regions"
echo "Candidates: $candidates"
bcftools stats "$candidates" | grep '^SN'
