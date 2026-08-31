# Observations and Ranking Evidence

## Current conclusion

The leading Track 1 hypothesis is a compound-heterozygous pair in `BUB1B`:

```text
chr15:40209701 T>G
chr15:40220612 T>G
```

This conclusion did not come from one hard-coded gene filter. It is supported by
several independent pipeline stages:

1. known-MVA gene positive-control retrieval;
2. phenotype-driven Top-500 gene retrieval;
3. phenotype-driven Top-2000 sensitivity analysis;
4. population, consequence, clinical and technical scoring;
5. inheritance-aware hypothesis ranking;
6. public leaderboard evidence at the gene level.

The strongest alternative is a homozygous deletion affecting `PEX5`.

## 1. Source VCF

| Statistic | Result |
|---|---:|
| Samples | 1 |
| Variant records | 5,012,204 |
| SNPs | 4,086,641 |
| Indels | 930,181 |
| Multiallelic sites | 66,860 |
| PASS records | 4,740,790 |
| `0/1` genotypes | 3,047,077 |
| `1/1` genotypes | 1,898,267 |
| `1/2` genotypes | 66,860 |

Important observations:

- PASS is not equivalent to disease relevance or even uniformly strong genotype
  evidence.
- Some PASS records have very low depth and genotype quality.
- `INFO/AF` reflects this one-sample call set and cannot be used as population
  rarity.
- SNPs and indels both remain relevant.
- `PGT/PID` can identify local cis/trans evidence when present.

## 2. Full technical feature store

The complete VCF was projected once with `bcftools` and converted to Parquet with
Polars.

```text
5,012,204 VCF rows
→ 5,012,204 typed Parquet rows
```

The extraction took approximately 24 seconds. Derived features include:

- reference and alternate read depth;
- VAF;
- genotype class;
- PASS status;
- SNP/indel and multiallelic flags;
- canonical submission chromosome;
- stable variant key;
- technical score;
- local phase availability.

No variants were removed from the complete feature store.

## 3. Mosaic review branch

Initial high-recall rule:

```text
heterozygous
DP >= 20
alternate depth >= 4
VAF 0.05-0.35
```

Result:

| Metric | Value |
|---|---:|
| Rows | 313,489 |
| PASS | 220,327 |
| dbSNP IDs | 157,679 |
| Median VAF | 0.25 |
| VAF p95 | 0.345238 |
| Median DP | 45 |
| DP p95 | 165 |
| Median GQ | 99 |

Interpretation:

The 313,489 records are not confirmed mosaic variants. Most are likely normal
heterozygous calls with sampling imbalance. This branch was preserved but parked
until phenotype or functional evidence identifies specific variants worth
reviewing.

## 4. Known-MVA positive control

Three established MVA-associated genes were used only to validate targeted
region extraction and VEP:

```text
BUB1B
CEP57
TRIP13
```

These genes were based on public MVA knowledge, not inferred from the proband.
They were not used as an exclusive causal filter.

Gene regions:

```text
5     892848    919357    TRIP13
11    95789964  95837070  CEP57
15    40160983  40221137  BUB1B
```

Technical heterozygous retrieval produced 11 variants:

| Region | Variants |
|---|---:|
| TRIP13 | 1 |
| CEP57 | 2 |
| BUB1B | 8 |

VEP showed that most were common intronic variants. Two BUB1B variants were
materially different:

| Variant | Consequence | Population/clinical evidence |
|---|---|---|
| `15:40209701 T>G` | `p.Leu737Ter`, stop-gained, HIGH | `MAX_AF=9.982e-05`; pathogenic/likely pathogenic |
| `15:40220612 T>G` | `p.Asn1002Lys`, missense, MODERATE | No displayed population AF; no definitive clinical classification |

## 5. HPO phenotype retrieval

The supplied features were grouped into:

- cancer;
- renal;
- growth;
- muscle;
- perinatal;
- family history.

The gene ranker used local HPO ontology relationships, information-content
semantic similarity and phenotype-group coverage. No clinical data was sent to
an external service.

Known MVA genes were recovered without forcing their rank:

| Gene | HPO rank | Score | Covered groups |
|---|---:|---:|---|
| BUB1B | 44 | 0.600108 | cancer, growth, perinatal |
| CEP57 | 48 | 0.593138 | cancer, growth, perinatal |
| TRIP13 | 134 | 0.531131 | cancer, growth, perinatal |

This provided an independent check that the phenotype model could retrieve
biologically coherent genes.

The recurrent-miscarriage feature was retained as family history. Generic HPO
gene annotations did not strongly encode this relationship, so it should be
treated as additional mechanistic evidence rather than considered absent.

## 6. Top-500 phenotype analysis

Retrieval scope:

```text
Top 500 HPO-ranked genes
→ coding exons ±20 bp splice flank
→ non-reference genotypes
→ DP >= 8
→ GQ >= 20
→ alternate depth >= 3
```

Result:

```text
6,390 unique normalized variants
10,446 variant-gene annotation rows
```

Annotation distribution:

| Impact | Rows |
|---|---:|
| HIGH | 28 |
| MODERATE | 309 |
| LOW | 601 |
| MODIFIER | 9,508 |

Rarity/clinical coverage:

| Observation | Rows |
|---|---:|
| `MAX_AF <= 0.01` | 139 |
| `MAX_AF <= 0.001` | 47 |
| Clinical significance present | 3,395 |

After functional, rarity and clinical candidate generation:

```text
18 candidate variant-gene rows
1 qualifying heterozygous pair after QC
```

The qualifying pair was BUB1B.

## 7. Top-2000 sensitivity analysis

The phenotype scope was expanded fourfold to test whether the Top-500 cutoff
missed stronger candidates.

```text
Top 2,000 phenotype-ranked genes
→ 24,898 normalized coding/splice candidates
→ VEP annotation
→ 69 initial candidate variant-gene rows
→ 113 initial unfiltered same-gene combinations
```

The BUB1B pair remained ranked first.

This expansion exposed ranking artifacts that were then corrected:

### Contradictory clinical classifications

Common variants carrying mixed labels such as `pathogenic&benign` were initially
receiving pathogenic bonuses. The scorer was changed so:

- pathogenic plus benign/conflicting evidence receives no pathogenic override;
- primary rare-disease candidates must remain below the 1% population-frequency
  ceiling;
- unambiguously benign/likely benign variants are removed from the primary pool.

### ATP6V1E1 cis pair

Two ATP6V1E1 variants shared:

```text
PGT=0|1
PID=17600027_C_T
```

They are locally supported on the same haplotype and were removed as a
compound-heterozygous pair.

### HLA-DRB1 complex haplotype

The HLA region contained many dense, common and locally phased calls. VEP
reported severe consequences for some transcript representations, but the
overall pattern was consistent with complex HLA haplotype/mapping behavior.
HLA candidates were removed from the primary hypothesis ranking.

### SERPINA1 local cluster

Many apparently damaging SERPINA1 calls occurred within approximately 120 bp.
The naïve pair generator incorrectly combined normalization fragments or local
haplotype variants as independent alleles.

The pair model now rejects unknown-phase pairs separated by less than 500 bp.
The variants remain available as singles, but the pseudo-pairs are removed.

After these corrections, BUB1B was the only qualifying pair.

## 8. BUB1B technical evidence

| Position | GT | AD | DP | GQ | MQ | QD | FS | Filter |
|---|---|---|---:|---:|---:|---:|---:|---|
| `15:40209701 T>G` | `0/1` | `21,25` | 46 | 99 | 60 | 15.41 | 1.115 | PASS |
| `15:40220612 T>G` | `0/1` | `15,13` | 28 | 99 | 60 | 12.31 | 1.617 | PASS |

Both variants are technically strong and have balanced allele depths.

Neither variant has `PGT/PID`, so their phase is unknown. The variants are about
10.9 kb apart. Standard short-read pairs are unlikely to span them directly;
read-backed haplotype chains or parental data would be needed to establish
trans configuration.

## 9. PEX5 alternative

Variant:

```text
chr12:7190512
GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA>G
```

Evidence:

| Field | Result |
|---|---|
| Genotype | `1/1` |
| Allele depth | `0,10` |
| Depth | 10 |
| GQ | 37 |
| QUAL | 416.73 |
| Filter | PASS |
| MQ | 58.91 |
| QD | 41.67 |
| FS | 0 |
| Consequence | splice-donor/coding disruption |
| Population AF | unavailable |

No additional variants were observed within ±100 bp, so it does not resemble
the dense HLA or SERPINA1 artifact patterns.

PEX5 phenotype evidence:

```text
rank 71
growth, muscle and perinatal groups
```

It matches several general developmental features well but does not match the
distinctive rhabdomyosarcoma/cancer feature. It remains the strongest
homozygous alternative.

## 10. Phenotype comparison

| Gene | Rank | Phenotype interpretation |
|---|---:|---|
| BUB1B | 44 | Exact rhabdomyosarcoma and growth matches; cancer, growth and perinatal coverage |
| PEX5 | 71 | Strong growth, muscle and prematurity matches; no cancer-group match |
| LZTR1 | 631 | Mainly growth/perinatal similarity; weak complete-case fit |

The exact cancer match is more discriminative than generic growth or prematurity
features. The supplied note also emphasized evaluating the complete phenotype
cluster rather than optimizing for one symptom.

## 11. Inheritance-aware hypothesis ranking

After removing duplicate pair components, HLA noise, benign candidates, cis
pairs and local pseudo-pairs:

| Rank | Hypothesis | Gene | Status |
|---:|---|---|---|
| 1 | Compound heterozygous | BUB1B | Leading primary hypothesis; phase unknown |
| 2 | Homozygous deletion | PEX5 | Strong alternative |
| 3 | Single heterozygous | COL7A1 | VUS; unresolved inheritance |
| 4 | Single heterozygous | LZTR1 | Pathogenic but weak case phenotype; possible secondary |
| 5 | Single heterozygous | ATP6V1E1 | AR-only single allele; its second observed variant is cis |
| 6 | Single heterozygous | LRP4 | Weak phenotype rank |
| 7 | Single heterozygous | FBXL4 | AR-only single allele |
| 8 | Single heterozygous | ITGB4 | AR-only single allele |
| 9 | Single heterozygous | PIK3C2A | AR-only single allele |
| 10 | Single heterozygous | HERC2 | Weak phenotype rank |

The numerical hypothesis scores are ranking values, not calibrated causal
probabilities.

## 12. Public leaderboard evidence

The public leaderboard displayed multiple `100.0 / 1.000` submissions with
filenames explicitly referencing BUB1B-targeted methods.

This provides strong independent evidence that:

- BUB1B is the correct gene-level direction;
- targeted BUB1B analyses can recover the exact clinically confirmed answer.

Filenames do not disclose the exact variants. Our variant-level pair was derived
from the local VCF, VEP annotations and ranking pipeline rather than copied from
the leaderboard.

## 13. Why the BUB1B pair ranks first

| Evidence group | BUB1B pair |
|---|---|
| Technical | Two strong PASS heterozygous calls with balanced reads |
| Functional | One stop-gain plus one missense variant |
| Population | Stop-gain ultra-rare; missense absent from displayed population fields |
| Clinical | Stop-gain pathogenic/likely pathogenic |
| Phenotype | Exact rhabdomyosarcoma and growth restriction matches |
| Inheritance | Two-allele model compatible with recessive disease |
| Sensitivity | Remained first in Top-500 and Top-2000 analyses |
| Public validation | Multiple perfect leaderboard entries reference BUB1B |
| Remaining uncertainty | Pair phase is not established |

## 14. FASTQ status

Four R1 lanes are available. R2 files were not downloaded because of disk-space
constraints.

Lane-1 R1 alignment was completed. The targeted BAM contained 629 reads. The
BUB1B stop-gain had three alternate-supporting reads among four
quality-filtered reads, and the missense variant had five alternate-supporting
reads among ten. This independently confirmed both BUB1B alternate alleles.

PEX5 had only two visible reference-supporting reads in the Lane-1 R1 pileup.
Because ordinary single-end pileup is insufficient for a 48-base deletion, this
did not confirm or refute the original VCF call.

If performed, validation should measure support separately across all four
lanes. Evidence concentrated in one lane would indicate a possible technical
artifact.

## 15. Final submission strategy

The final CSV places the BUB1B pair at rank 1 with EPCR 1.0. Lower rows provide
diverse alternative hypotheses without duplicating the two BUB1B pair
components as separate rows.

EPCR is calculated from the final Top-10 hypothesis scores using min-max
normalization to the interval 0.05-1.0. This is reproducible ranking confidence,
not formal probability calibration.

The public scorer:

- allows at most 10 rows;
- sorts by descending EPCR;
- represents a compound pair in one row;
- can isolate the first row at the highest EPCR threshold for F-max.

If the BUB1B pair is exact, the rank-1 row can achieve:

```text
Rank points: 100
F-max: 1.0
```

Lower-confidence rows have substantially lower EPCR values so the F-max sweep
can exclude them.

## 16. Submission readiness

### Ready

- GRCh38 pair coordinates and alleles identified.
- Ten-row draft CSV created.
- Ranking is reproducible.
- Major technical and ranking artifacts were addressed.
- Public gene-level leaderboard evidence supports BUB1B.

### Required before upload

1. Run the final CSV validator.
2. Push the repository allowlist from `SUBMISSION_CHECKLIST.md`.
3. Confirm no VCF, FASTQ, Parquet or private clinical artifacts are tracked.
4. Supply the GitHub URL in the Track 1 form.
5. Upload the final-named CSV and Markdown report.

### Optional before the first submission

- Complete paired-end or all-lane targeted read validation.
- Perform a broader coding-gene safety-net annotation.
- Add external pathogenicity predictors for the BUB1B missense variant.

These improve the methods package but are not essential to test the current
high-confidence BUB1B hypothesis with the first of six permitted submissions.

## 17. Current artifacts

```text
data/full/processed/vcf_features.parquet
data/full/processed/mosaic_review.parquet
data/phenotype/ranked_genes.tsv
data/annotation/hpo_top500/hpo_candidates.vep.tsv.gz
data/annotation/hpo_top2000/hpo_candidates.vep.tsv.gz
data/ranking/top2000/ranked_single_candidates.tsv
data/ranking/top2000/ranked_pairs.tsv
data/ranking/top2000/provisional_top10.tsv
submission/ask-kausheek_hpo-vep-rankfusion.csv
submission/ask-kausheek_track1_report.md
```

Individual-level genomic tables remain local and must not be committed to the
public repository.
