# C. Pipeline design

```text
VCF.gz
  -> one-pass bcftools projection
  -> typed Parquet feature table
  -> normalization and local annotation joins
  -> four retrieval pools
  -> phenotype/gene scoring
  -> same-gene compound-het pair generation
  -> interpretable component and hypothesis scores
  -> top 10 pair/single rows
  -> targeted R1 validation for leading alleles
  -> submission CSV
```

## Stages

| Stage | Transformation | Search-space effect | Persisted output |
|---|---|---|---|
| Ingestion | Project required VCF INFO/FORMAT columns once | No filtering | `vcf_features.parquet` |
| Normalization | Split multiallelics and left-normalize before external joins | Stable variant keys | normalized VCF/Parquet |
| Annotation | Local joins for consequence, population, ClinVar | Adds ranking signals | `annotated_variants.parquet` |
| Retrieval A | Rare, technically credible, functional variants | Millions to thousands | `pool_rare.parquet` |
| Retrieval B | Known pathogenic/likely pathogenic variants | Independent safety net | `pool_clinvar.parquet` |
| Retrieval C | Low-VAF/mosaic candidates | Separate nonstandard-genotype branch | `mosaic_review.parquet` |
| Retrieval D | CNV/SV/chromosome evidence | Separate event representation | event table |
| Pairing | Same-gene heterozygous combinations with inheritance support | Thousands to manageable pairs | `candidate_pairs.parquet` |
| Ranking | Weighted normalized features plus inheritance-aware hypothesis adjustments | Top tens | ranked table |
| Validation | Read-level checks only for finalists | Removes artifacts | evidence table |

## Annotation decision gate

Do not wait indefinitely for the full VEP cache. Benchmark the fastest locally permitted route on chromosome 22:

1. Measure wall time and annotation coverage.
2. Confirm stable GRCh38 variant-key matching.
3. Continue only if projected whole-genome runtime is acceptable.

Prefer a local database join when it supplies the minimum feature set faster than transcript-heavy annotation.

## Progressive VEP scope

1. Positive control: literature-derived MVA genes from
   `configs/known_mva_seed_genes.txt`; verifies extraction and annotation.
2. Primary retrieval: top 500 genes from case HPO terms.
3. Sensitivity retrieval: top 2,000 phenotype-ranked genes.
4. Whole-genome annotation only if earlier stages fail and measured throughput
   is acceptable.

The three-gene positive control is not a proband-derived filter and cannot by
itself establish causality. The primary analysis starts with phenotype-ranked
genes and retains explicit expansion paths.

## Branch combination

The pair branch receives priority because the public scorer explicitly describes
a compound-heterozygous production answer. It must still compete against
homozygous, single-variant, mosaic and CNV/SV hypotheses using data-derived
evidence. No pair should rank first solely because of scorer metadata.

## Implemented sensitivity checks

- Top-500 versus Top-2000 phenotype retrieval
- Clinical-prior versus equal HPO feature weights
- Missing population AF changed from favorable to neutral
- Cis pair exclusion using `PGT/PID`
- Unknown-phase local-cluster suppression
- HLA dense-haplotype exclusion
- One single-variant fallback per gene

BUB1B remained the leading pair across these checks.
