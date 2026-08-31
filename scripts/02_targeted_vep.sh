#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 6 || $# -gt 7 ]]; then
  cat >&2 <<'EOF'
Usage:
  02_targeted_vep.sh INPUT.vcf.gz REFERENCE.fa GENES.txt GENES.gtf.gz \
    VEP_CACHE_DIR OUTPUT_DIR [THREADS]
EOF
  exit 2
fi

input=$1
reference=$2
genes_file=$3
gtf=$4
cache_dir=$5
output_dir=$6
threads=${7:-8}

for command in bcftools bgzip tabix vep awk sort gzip comm grep sed cut; do
  command -v "$command" >/dev/null || {
    echo "Missing required command: $command" >&2
    exit 1
  }
done

for file in "$input" "$reference" "$genes_file" "$gtf"; do
  [[ -f "$file" ]] || {
    echo "Required file not found: $file" >&2
    exit 1
  }
done

mkdir -p "$output_dir"
regions="$output_dir/gene_regions.bed"
subset="$output_dir/targeted_candidates.vcf.gz"
annotated="$output_dir/targeted_candidates.vep.vcf"

# Convert selected GTF gene records from 1-based closed coordinates to BED.
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

missing=$(
  comm -23 \
    <(grep -v '^#' "$genes_file" | sed '/^$/d' | sort -u) \
    <(cut -f4 "$regions" | sort -u)
)
if [[ -n "$missing" ]]; then
  echo "Genes missing from GTF:" >&2
  echo "$missing" >&2
  exit 1
fi

# Region restriction happens before normalization, so only a tiny fraction of
# the genome reaches VEP. Borderline FILTER values remain eligible.
bcftools view -R "$regions" "$input" -Ou |
  bcftools norm -f "$reference" -m -any --keep-sum AD -Ou |
  bcftools view \
    -i 'GT="het" && FMT/DP>=8 && FMT/GQ>=20 && FMT/AD[0:1]>=3' \
    -Oz -o "$subset"
tabix -f -p vcf "$subset"

vep \
  --input_file "$subset" \
  --output_file "$annotated" \
  --format vcf \
  --vcf \
  --offline \
  --cache \
  --dir_cache "$cache_dir" \
  --assembly GRCh38 \
  --fasta "$reference" \
  --fork "$threads" \
  --buffer_size 5000 \
  --symbol \
  --canonical \
  --mane \
  --biotype \
  --variant_class \
  --hgvs \
  --numbers \
  --total_length \
  --af \
  --af_1kg \
  --af_gnomade \
  --af_gnomadg \
  --max_af \
  --check_existing \
  --pick_allele_gene \
  --stats_file "$output_dir/targeted_candidates.vep.summary.html" \
  --force_overwrite

bgzip -f "$annotated"
tabix -f -p vcf "${annotated}.gz"

echo "Regions:   $regions"
echo "Candidates: $subset"
echo "Annotated:  ${annotated}.gz"
