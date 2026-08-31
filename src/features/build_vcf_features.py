from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl
import yaml


STRING_COLUMNS = {
    "chrom": pl.String,
    "id": pl.String,
    "ref": pl.String,
    "alt": pl.String,
    "filter": pl.String,
    "gt": pl.String,
    "ad": pl.String,
    "pl": pl.String,
    "pgt": pl.String,
    "pid": pl.String,
}


def load_thresholds(config_path: Path) -> dict[str, float]:
    with config_path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)["retrieval"]


def build_features(input_path: Path) -> pl.LazyFrame:
    frame = pl.scan_csv(
        input_path,
        separator="\t",
        null_values=".",
        schema_overrides=STRING_COLUMNS,
        infer_schema_length=50_000,
    )

    ref_depth = pl.col("ad").str.split(",").list.get(0, null_on_oob=True).cast(
        pl.Int32, strict=False
    )
    sample_depth = pl.col("sample_dp").cast(pl.Int32, strict=False)
    alt_depth = (sample_depth - ref_depth).clip(lower_bound=0)
    submission_chrom = (
        pl.when(pl.col("chrom").str.starts_with("chr"))
        .then(pl.col("chrom"))
        .otherwise(pl.lit("chr") + pl.col("chrom"))
    )

    return frame.with_columns(
        pl.col("pos").cast(pl.Int64, strict=False),
        pl.col("qual").cast(pl.Float32, strict=False),
        pl.col("gq").cast(pl.Int16, strict=False),
        sample_depth.alias("sample_dp"),
        ref_depth.alias("ref_depth"),
        alt_depth.alias("alt_depth"),
        submission_chrom.alias("submission_chrom"),
        pl.concat_str(
            [submission_chrom, pl.col("pos"), pl.col("ref"), pl.col("alt")],
            separator=":",
        ).alias("variant_key"),
        pl.col("id").str.starts_with("rs").fill_null(False).alias("is_dbsnp"),
        (pl.col("filter") == "PASS").alias("is_pass"),
        pl.col("alt").str.contains(",").alias("is_multiallelic"),
        ((pl.col("ref").str.len_chars() == 1) & (pl.col("alt").str.len_chars() == 1)).alias(
            "is_snp"
        ),
        pl.col("gt").is_in(["0/1", "1/0", "0|1", "1|0"]).alias("is_het_ref_alt"),
        pl.col("gt").is_in(["1/1", "1|1"]).alias("is_hom_alt"),
        (pl.col("pid").is_not_null() & pl.col("pgt").is_not_null()).alias(
            "has_local_phase"
        ),
    ).with_columns(
        pl.when(pl.col("sample_dp") > 0)
        .then(pl.col("alt_depth") / pl.col("sample_dp"))
        .otherwise(None)
        .cast(pl.Float32)
        .alias("vaf"),
        (
            pl.col("qual").fill_null(0).clip(0, 200) / 200 * 0.25
            + pl.col("gq").fill_null(0).clip(0, 99) / 99 * 0.35
            + pl.col("sample_dp").fill_null(0).clip(0, 30) / 30 * 0.25
            + pl.when(pl.col("is_pass")).then(0.15).otherwise(0.0)
        )
        .cast(pl.Float32)
        .alias("technical_score"),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/proband.yaml"),
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    thresholds = load_thresholds(args.config)
    features = build_features(args.input)

    features.sink_parquet(
        args.output_dir / "vcf_features.parquet",
        compression="zstd",
        statistics=True,
    )

    mosaic = features.filter(
        pl.col("is_het_ref_alt")
        & (pl.col("sample_dp") >= thresholds["mosaic_minimum_depth"])
        & (pl.col("alt_depth") >= thresholds["mosaic_minimum_alt_depth"])
        & pl.col("vaf").is_between(
            thresholds["mosaic_vaf_min"],
            thresholds["mosaic_vaf_max"],
            closed="both",
        )
    )
    mosaic.sink_parquet(
        args.output_dir / "mosaic_review.parquet",
        compression="zstd",
        statistics=True,
    )


if __name__ == "__main__":
    main()
