#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 INPUT.vcf.gz OUTPUT_DIR" >&2
  exit 2
fi

input=$1
output_dir=$2
mkdir -p "$output_dir"

for command in bcftools bgzip; do
  command -v "$command" >/dev/null || {
    echo "Missing required command: $command" >&2
    exit 1
  }
done

[[ -f "$input" ]] || {
  echo "VCF not found: $input" >&2
  exit 1
}
[[ -f "${input}.tbi" || -f "${input}.csi" ]] || {
  echo "VCF index not found beside: $input" >&2
  exit 1
}

stats_file="$output_dir/bcftools.stats.txt"
features_file="$output_dir/vcf_features.tsv.gz"

bcftools stats "$input" >"$stats_file"

{
  printf '%s\n' $'chrom\tpos\tid\tref\talt\tqual\tfilter\tac\taf\tan\tinfo_dp\tfs\tmq\tqd\tmq_rank_sum\tread_pos_rank_sum\tbase_q_rank_sum\tclipping_rank_sum\tsor\texcess_het\tinbreeding_coeff\tgt\tad\tsample_dp\tgq\tpl\tpgt\tpid'
  bcftools query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%INFO/DP\t%INFO/FS\t%INFO/MQ\t%INFO/QD\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\t%INFO/BaseQRankSum\t%INFO/ClippingRankSum\t%INFO/SOR\t%INFO/ExcessHet\t%INFO/InbreedingCoeff[\t%GT\t%AD\t%DP\t%GQ\t%PL\t%PGT\t%PID]\n' \
    "$input"
} | bgzip -c >"$features_file"

printf 'Created %s\nCreated %s\n' "$stats_file" "$features_file"
