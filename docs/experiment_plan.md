# E. Experiment plan

## Offline benchmark

Build a public, non-case-specific benchmark of known pathogenic compound-heterozygous cases plus matched benign/background variants. Split by gene and disease to prevent variant/gene leakage.

Track:

| Metric | Reason |
|---|---|
| Top-1/5/10 pair recall | Matches leaderboard objective |
| MRR and NDCG@10 | Measures ranking quality |
| Candidate count and annotation coverage | Detects over-filtering |
| Runtime and peak memory | Keeps iteration practical |

## Completed case sensitivity experiments

| ID | Change |
|---|---|
| E00 | Complete technical feature extraction |
| E01 | Known-MVA positive-control retrieval |
| E02 | HPO Top-500 coding/splice retrieval |
| E03 | HPO Top-2000 sensitivity expansion |
| E04 | Clinical-prior HPO weights |
| E05 | Equal HPO weights |
| E06 | Missing population AF favorable versus neutral |
| E07 | Pair QC: cis, local clusters and HLA artifacts |
| E08 | Lane-1 R1 confirmation of BUB1B alleles |

BUB1B remained the leading pair in the Top-500, Top-2000 and equal-HPO-weight
experiments. PEX5 remained the leading homozygous alternative.

## Submission budget

Do not spend a submission on technical-only rankings. Suggested gates:

1. Submission 1: best validated baseline ensemble.
2. Submission 2: phenotype-heavy ablation only if top pairs materially differ.
3. Submission 3: pathogenicity-heavy ablation.
4. Submission 4: mosaic-aware reranking if read evidence supports it.
5. Submission 5: model/rank-fusion ensemble after error analysis.
6. Submission 6: reserve for the final validated run.

Record each run in `experiments/experiments.csv` with data hash, annotation versions, thresholds, features, model/weights, candidate counts, runtime, and leaderboard feedback.

## Remaining evaluation work

The private causal label prevents offline Top-K measurement for this proband.
The immediate evaluation is the live Submission-1 score:

1. Upload the validated CSV, methods report and repository URL.
2. Record rank points, F-max and full/partial/no-match feedback.
3. If the result is perfect, freeze the Track 1 prediction.
4. If not, use the remaining submissions for broader coding, mosaic/CNV or
   additional read-validation experiments.

Future method development can use public solved rare-disease cases with
gene/disease-held-out splits for Top-K, MRR and NDCG evaluation.
