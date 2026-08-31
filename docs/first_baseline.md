# D. First baseline

The baseline is an interpretable retrieval-and-ranking system, not a trained model.

## Candidate retrieval

1. Normalize variants and split multiallelic records.
2. Keep a broad main pool with adequate genotype support; do not require PASS for every branch.
3. Retrieve rare functional variants using local population and consequence annotations.
4. Retrieve ClinVar pathogenic/likely pathogenic variants independently.
5. Rank HPO-relevant genes.
6. Generate pairs only when both variants map to the same gene and are compatible with a compound-heterozygous hypothesis.

The first VEP run over three known MVA genes is a pipeline positive control. The
primary case analysis must repeat the same workflow over phenotype-ranked genes,
regardless of whether the positive-control set contains an interesting pair.

## Implemented component scores

Single candidates combine normalized:

- phenotype relevance;
- functional consequence;
- population rarity;
- clinical evidence;
- technical evidence.

Pair candidates additionally require both alleles to satisfy rarity/pathogenicity
criteria and include combined consequence, technical and phase evidence.

The weights are heuristic experiment parameters, not learned biological truth.
The leading result was therefore tested with equal HPO feature weights and a
neutral missing-AF score. BUB1B remained first.

For variants sharing a `PID`, use `PGT` to determine local haplotype assignment. Do not infer trans configuration merely because two heterozygous calls occur in the same gene; unphased finalists require read-level or parental evidence.

## Technical feature layer

`scripts/01_extract_vcf_features.sh` makes one reusable compressed projection. `src/features/build_vcf_features.py` converts it to Parquet and creates `mosaic_review.parquet`.

Technical fields alone do not emit final causal candidates; they are joined with
VEP, phenotype and inheritance evidence.

## Submission-aware ranking

Official rank tiers are:

| Rank | Points |
|---:|---:|
| 1 | 100 |
| 2-3 | 50 |
| 4-5 | 25 |
| 6-10 | 10 |

The exact pair must occupy one row for full credit. EPCR values control sorting
and F-max thresholds. Generate them reproducibly from the final hypothesis
scores while documenting that they are not calibrated probabilities. Put the
strongest pair first and avoid filling rows with duplicate hypotheses.

The scorer currently includes all rows in F-max regardless of `finding_type`. Keep secondary findings below primary candidates with lower EPCR so a threshold can exclude them.
