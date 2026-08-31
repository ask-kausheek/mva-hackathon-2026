#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 REFERENCE.fa VEP_CACHE_DIR OUTPUT_DIR [THREADS]" >&2
  exit 2
fi

reference=$1
cache_dir=$2
output_dir=$3
threads=${4:-8}

work_mount=$(realpath "$output_dir")
cache_mount=$(realpath "$cache_dir")
reference_mount=$(realpath "$(dirname "$reference")")
reference_name=$(basename "$reference")
input="$output_dir/targeted_candidates.vcf.gz"

[[ -f "$input" ]] || {
  echo "Run 02_prepare_targeted_candidates.sh first; missing $input" >&2
  exit 1
}

docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "$work_mount:/work" \
  -v "$cache_mount:/opt/vep/.vep:ro" \
  -v "$reference_mount:/reference:ro" \
  ensemblorg/ensembl-vep \
  vep \
  --input_file /work/targeted_candidates.vcf.gz \
  --output_file /work/targeted_candidates.vep.vcf \
  --format vcf \
  --vcf \
  --offline \
  --cache \
  --dir_cache /opt/vep/.vep \
  --assembly GRCh38 \
  --fasta "/reference/$reference_name" \
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
  --no_stats \
  --force_overwrite

bgzip -f "$output_dir/targeted_candidates.vep.vcf"
tabix -f -p vcf "$output_dir/targeted_candidates.vep.vcf.gz"
echo "Created $output_dir/targeted_candidates.vep.vcf.gz"
