# Observations and Ranking Evidence

## Result

The leading hypothesis is a compound-heterozygous BUB1B pair:

| Variant | Genotype | Consequence | Key evidence |
|---|---|---|---|
| `chr15:40209701 T>G` | `0/1` | `p.Leu737Ter`, stop-gained | Ultra-rare; pathogenic/likely pathogenic |
| `chr15:40220612 T>G` | `0/1` | `p.Asn1002Lys`, missense | No displayed population frequency |

The strongest competing model is a homozygous deletion in PEX5 at
`chr12:7190512`.

## Source-data profile

| Statistic | Value |
|---|---:|
| Variant records | 5,012,204 |
| SNPs | 4,086,641 |
| Indels | 930,181 |
| Multiallelic sites | 66,860 |
| PASS records | 4,740,790 |
| Heterozygous genotypes | 3,047,077 |
| Homozygous alternate genotypes | 1,898,267 |

The VCF was projected once and converted into a typed Parquet feature store
without removing variants. This avoided repeated scans during candidate
retrieval and sensitivity analysis.

`INFO/AF` describes the one-sample call set and was not used as population
frequency. Population rarity came from annotation resources.

## Technical feature store

Derived features included:

- normalized variant identity;
- genotype class;
- reference and alternate depth;
- sample depth and genotype quality;
- VAF;
- filter and variant-type flags;
- technical-quality score;
- local phase fields.

The initial mosaic-review branch selected heterozygous calls with DP at least
20, alternate depth at least 4 and VAF from 0.05 through 0.35.

| Mosaic-review metric | Value |
|---|---:|
| Candidate rows | 313,489 |
| PASS rows | 220,327 |
| Median VAF | 0.25 |
| Median DP | 45 |
| Median GQ | 99 |

The distribution was consistent with a high-recall pool dominated by ordinary
heterozygous imbalance. It was retained as a review branch rather than promoted
directly into the final ranking.

## Phenotype retrieval

The local HPO ranker used semantic similarity, information content and coverage
across phenotype groups.

| Gene | HPO rank | Score | Covered groups |
|---|---:|---:|---|
| BUB1B | 44 | 0.600108 | Cancer, growth, perinatal |
| CEP57 | 48 | 0.593138 | Cancer, growth, perinatal |
| PEX5 | 71 | — | Growth, muscle, perinatal |
| TRIP13 | 134 | 0.531131 | Cancer, growth, perinatal |

BUB1B was recovered by phenotype ranking without forcing it to the top. PEX5
matched several general developmental features but lacked the stronger
cancer-group match.

## Retrieval experiments

| Experiment | Normalized candidates | Candidate rows after scoring | Initial pair combinations | Leading result |
|---|---:|---:|---:|---|
| Top-500 genes | 6,390 | 18 | 1 qualifying pair | BUB1B pair |
| Top-2000 genes | 24,898 | 69 | 113 before pair QC | BUB1B pair |

Expanding retrieval fourfold did not displace the leading pair. The larger pool
was useful for exposing ranking and pair-generation artifacts.

## Ranking corrections

### Clinical-label conflicts

Variants carrying both pathogenic and benign evidence initially received an
incorrect pathogenic bonus. Conflicting classifications now receive no
pathogenic override, and unambiguously benign variants are excluded from the
primary pool.

### Missing population frequency

Missing frequency initially behaved like positive rarity evidence. It is now
assigned a neutral score so absence of annotation does not imply rarity.

### Known-cis pair

Two ATP6V1E1 variants shared `PGT=0|1` and
`PID=17600027_C_T`. They were excluded as a compound-heterozygous pair because
the local phase fields place them on the same haplotype.

### Local pseudo-pairs

Closely spaced SERPINA1 calls produced pair combinations that were more
consistent with local haplotype or normalization artifacts than independent
alleles. Unknown-phase pairs within 500 bases are now excluded.

### Complex-region calls

Dense HLA-DRB1 calls were removed from primary ranking because transcript-level
severe consequences occurred within a common, complex haplotype pattern.

After these controls, BUB1B remained the only qualifying heterozygous pair.

## BUB1B evidence

| Field | `chr15:40209701 T>G` | `chr15:40220612 T>G` |
|---|---:|---:|
| Genotype | `0/1` | `0/1` |
| Allele depth | `21,25` | `15,13` |
| Depth | 46 | 28 |
| GQ | 99 | 99 |
| MQ | 60 | 60 |
| QD | 15.41 | 12.31 |
| FS | 1.115 | 1.617 |
| Filter | PASS | PASS |

Both calls have balanced allele depths and strong genotype and mapping quality.
They are approximately 10.9 kb apart and have no shared `PGT/PID`; phase is
therefore unresolved.

Evidence convergence:

| Evidence dimension | Observation |
|---|---|
| Technical | Two strong heterozygous PASS calls |
| Functional | Stop-gain plus missense consequence |
| Population | Ultra-rare stop-gain; missing displayed AF for missense |
| Clinical | Stop-gain classified pathogenic/likely pathogenic |
| Phenotype | Strong cancer, growth and perinatal coverage |
| Inheritance | Compatible two-allele recessive model |
| Sensitivity | First under both retrieval scopes and HPO-weight settings |

## PEX5 alternative

The PEX5 candidate is a homozygous 48-base deletion:

```text
chr12:7190512
GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA>G
```

| Field | Value |
|---|---|
| Genotype | `1/1` |
| Allele depth | `0,10` |
| Depth | 10 |
| GQ | 37 |
| QUAL | 416.73 |
| MQ | 58.91 |
| QD | 41.67 |
| FS | 0 |
| Filter | PASS |
| Consequence | Splice/coding disruption |
| Population frequency | Unavailable |

It remains the strongest homozygous alternative. Its phenotype evidence is
broader and less distinctive than the BUB1B match.

## Sensitivity analysis

The ranking was rerun with all HPO feature weights set to 1.0 and missing
population frequency treated neutrally. BUB1B remained first and PEX5 remained
second.

This result reduces dependence on:

- the original clinical weighting heuristic;
- the Top-500 retrieval cutoff;
- favorable treatment of missing annotations;
- unfiltered pair combinations.

## Targeted read validation

Lane-1 R1 alignment produced a targeted BAM containing 629 reads.

| Candidate | Quality-filtered support |
|---|---:|
| BUB1B stop-gain | 3 of 4 reads supported ALT |
| BUB1B missense | 5 of 10 reads supported ALT |

Both submitted BUB1B alleles were independently observed in Lane 1.

The PEX5 locus showed two reference-supporting reads in ordinary pileup. This
does not confirm or reject a 48-base deletion because the available single-end
pileup is not an adequate indel-validation method.

## Final ranked hypotheses

| Rank | Gene | Model | EPCR |
|---:|---|---|---:|
| 1 | BUB1B | Compound heterozygous | 1.000 |
| 2 | PEX5 | Homozygous deletion | 0.650 |
| 3 | LZTR1 | Single heterozygous; secondary | 0.428 |
| 4 | COL7A1 | Single heterozygous VUS | 0.330 |
| 5 | GNRHR | Single heterozygous allele in AR gene | 0.200 |
| 6 | ATP6V1E1 | Single heterozygous allele in AR gene | 0.145 |
| 7 | FAT4 | Single heterozygous allele in AR gene | 0.118 |
| 8 | FBXL4 | Single heterozygous allele in AR gene | 0.096 |
| 9 | SERPINA1 | Single heterozygous allele in AR gene | 0.076 |
| 10 | PIK3C2A | Single heterozygous allele in AR gene | 0.050 |

EPCR values are normalized ranking confidence, not calibrated causal
probabilities.

## Remaining uncertainty

- BUB1B phase has not been established.
- Parental segregation data are unavailable.
- Validation used one R1 lane rather than paired reads across all lanes.
- The PEX5 deletion requires an indel-aware validation method.
- CNV and structural-variant hypotheses were not implemented as a complete
  genome-wide branch.
