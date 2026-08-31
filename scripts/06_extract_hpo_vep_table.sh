#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 ANNOTATED.vcf.gz OUTPUT_DIR" >&2
  exit 2
fi

input=$1
output_dir=$2
output="$output_dir/hpo_candidates.vep.tsv.gz"

mkdir -p "$output_dir"

{
  printf '%s\n' $'chrom\tpos\tref\talt\tfilter\tgt\tad\tdp\tgq\tpgt\tpid\tgene\tconsequence\timpact\thgvsc\thgvsp\trsid\tgnomade_af\tgnomadg_af\tmax_af\tclin_sig\tcanonical\tmane'

  bcftools +split-vep "$input" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\t[%AD]\t[%DP]\t[%GQ]\t[%PGT]\t[%PID]\t%SYMBOL\t%Consequence\t%IMPACT\t%HGVSc\t%HGVSp\t%Existing_variation\t%gnomADe_AF\t%gnomADg_AF\t%MAX_AF\t%CLIN_SIG\t%CANONICAL\t%MANE_SELECT\n' \
    -d
} | bgzip -c >"$output"

echo "Created $output"

python - "$output" <<'PY'
import sys
import polars as pl

path = sys.argv[1]
df = pl.read_csv(
    path,
    separator="\t",
    null_values=".",
    infer_schema_length=50_000,
)

unique_variants = df.select(
    pl.concat_str(["chrom", "pos", "ref", "alt"], separator=":")
      .n_unique()
      .alias("unique_variants"),
    pl.len().alias("annotation_rows"),
)
print("\nCounts")
print(unique_variants)

for column in ("gt", "impact"):
    print(f"\n{column}")
    print(df.group_by(column).len().sort("len", descending=True))

print("\nTop consequences")
print(
    df.group_by("consequence")
      .len()
      .sort("len", descending=True)
      .head(20)
)

print("\nPopulation-frequency availability")
print(
    df.select(
        pl.col("max_af").is_not_null().sum().alias("with_max_af"),
        (pl.col("max_af").fill_null(1.0) <= 0.01).sum().alias("max_af_le_1pct"),
        (pl.col("max_af").fill_null(1.0) <= 0.001).sum().alias("max_af_le_0_1pct"),
        pl.col("clin_sig").is_not_null().sum().alias("with_clin_sig"),
    )
)
PY
