# Track 1 Methods Report

**Participant:** ask-kausheek  
**Challenge:** Rare Disease, Real Kid: MVA Hackathon 2026  
**Track:** Track 1 — Ranked Variant Prediction  
**Reference assembly:** GRCh38  
**Proband:** PROBAND01

## Executive summary

I treated variant identification as a progressive retrieval and ranking problem
rather than applying one irreversible variant filter. A complete technical
feature store was created from 5,012,204 VCF records. Candidate retrieval then
used locally computed HPO semantic similarity, coding/splice-region selection,
VEP consequence and population annotations, technical evidence, inheritance
compatibility, and pair-level ranking.

The leading prediction is a two-variant `BUB1B` hypothesis:

- `chr15:40209701 T>G`, stop-gained, `p.Leu737Ter`
- `chr15:40220612 T>G`, missense, `p.Asn1002Lys`

Both variants are heterozygous PASS calls with balanced allele depths and GQ 99.
The stop-gain is ultra-rare and classified pathogenic/likely pathogenic. The
missense variant is absent from the displayed population-frequency fields.
Their phase is unresolved.

The strongest alternative is a homozygous splice/coding deletion in `PEX5` at
`chr12:7190512`.

## Input data

- One-sample GRCh38 VCF produced by Sentieon/GATK-compatible tooling
- 5,012,204 variant records
- 4,086,641 SNPs
- 930,181 indels
- 66,860 multiallelic sites
- 4,740,790 PASS records
- Four sequencing lanes; R1 FASTQs were locally available during the analysis
- Eight supplied HPO terms describing cancer, renal, growth, muscle, perinatal
  and family-history features

Raw genomic files and full individual-level derived tables remained local and
were excluded from the repository.

## Pipeline

```text
VCF
→ one-pass technical feature extraction
→ Parquet feature store
→ local HPO semantic gene ranking
→ coding/splice candidate retrieval
→ VEP annotation
→ rarity/function/clinical scoring
→ inheritance-aware single and pair hypotheses
→ Top-500/Top-2000 sensitivity analysis
→ ranked Track 1 CSV
```

## Technical feature engineering

`bcftools` streamed the VCF once and projected the required INFO and FORMAT
fields. Polars converted the projection into typed Parquet and derived:

- genotype class;
- reference and alternate read depth;
- VAF;
- depth and genotype quality;
- PASS status;
- SNP/indel and multiallelic flags;
- canonical submission chromosome;
- stable variant key;
- local phase availability from `PGT/PID`.

The complete Parquet retained all 5,012,204 variants.

## Mosaic review branch

A high-recall mosaic review pool selected heterozygous calls with:

- DP at least 20;
- alternate depth at least 4;
- VAF from 0.05 through 0.35.

This produced 313,489 records. The median VAF was 0.25, median GQ was 99,
220,327 records were PASS, and 157,679 had dbSNP identifiers. This indicated
that most of the pool represented ordinary heterozygous allele imbalance rather
than strong mosaic evidence. The pool was retained as a parallel safety net but
was not used as the primary ranking universe.

## Phenotype-driven gene retrieval

The HPO terms were processed locally using the public HPO ontology and
gene-to-phenotype annotations. The baseline ranker:

1. propagated annotations through HPO ancestor relationships;
2. computed information-content semantic similarity;
3. applied feature weights;
4. rewarded coverage across cancer, renal, growth, muscle, perinatal and
   family-history groups.

The family history of recurrent pregnancy loss was retained as a weighted
phenotype input.

Known MVA genes were recovered without forcing their positions:

| Gene | HPO rank |
|---|---:|
| BUB1B | 44 |
| CEP57 | 48 |
| TRIP13 | 134 |

## Candidate retrieval and VEP annotation

The primary retrieval selected exons plus 20-base splice flanks for
phenotype-ranked genes. It retained technically credible non-reference genotypes
with DP at least 8, GQ at least 20, and alternate depth at least 3.

### Top-500 experiment

- 6,390 unique normalized candidates
- 10,446 variant-gene annotation rows
- 28 HIGH-impact rows
- 309 MODERATE-impact rows
- 139 rows with maximum population AF at most 1%
- 47 rows with maximum population AF at most 0.1%

### Top-2000 sensitivity experiment

- 24,898 normalized coding/splice candidates
- 69 initial rare functional variant-gene candidates
- 113 initial same-gene combinations before pair QC

VEP release 116 supplied gene, transcript, consequence, impact, HGVS,
canonical/MANE status, existing identifiers, population frequencies and
clinical significance where present.

## Ranking

Single-variant component scores combined:

- phenotype relevance;
- functional consequence;
- population rarity;
- clinical evidence;
- technical evidence.

Pair scoring additionally considered:

- two qualifying heterozygous alleles in one gene;
- combined rarity and functional severity;
- recessive inheritance compatibility;
- phase evidence where available.

The following artifacts were explicitly controlled:

- cis-supported pairs were excluded;
- unknown-phase variants less than 500 bases apart were not treated as
  independent compound-heterozygous alleles;
- unambiguously benign variants were excluded;
- conflicting benign/pathogenic labels did not receive a pathogenic override;
- dense HLA haplotype calls were removed from primary hypothesis ranking;
- pair components were not duplicated in the final hypothesis Top-10.

The numerical scores were used for ordering and were not treated as calibrated
causal probabilities. Missing population frequency was assigned a neutral
rarity score rather than being treated as evidence of rarity.

Final EPCR values were generated reproducibly from the selected Top-10
hypothesis scores:

```text
EPCR = 0.05 + 0.95 ×
       (score - minimum Top-10 score) /
       (maximum Top-10 score - minimum Top-10 score)
```

This monotonic transformation maps the leading hypothesis to 1.0 and the lowest
submitted hypothesis to 0.05. It preserves ranking and threshold separation but
is not claimed as formal probability calibration.

## Leading BUB1B hypothesis

| Feature | `chr15:40209701 T>G` | `chr15:40220612 T>G` |
|---|---|---|
| Genotype | 0/1 | 0/1 |
| Allele depth | 21,25 | 15,13 |
| Depth | 46 | 28 |
| GQ | 99 | 99 |
| MQ | 60 | 60 |
| QD | 15.41 | 12.31 |
| Filter | PASS | PASS |
| Consequence | stop-gained | missense |
| Protein | p.Leu737Ter | p.Asn1002Lys |
| Population evidence | MAX_AF 9.982e-05 | no displayed AF |
| Clinical evidence | pathogenic/likely pathogenic | unresolved |

The pair ranked first in both the Top-500 and Top-2000 phenotype experiments.
It was the only pair remaining after rarity, phase, local-cluster and technical
QC.

The two variants are approximately 10.9 kb apart and have no `PGT/PID`, so phase
could not be established from the VCF.

## PEX5 alternative

The strongest alternative is:

```text
chr12:7190512
GGCCTCTGAGGCAGTGAGTGTTCTTGAGGTGGAAAGCCCAGGTGCA>G
```

It is a homozygous PASS deletion with 10/10 alternate reads, GQ 37, MQ 58.91
and a predicted splice/coding consequence. It is technically credible and
compatible with a recessive model. Its phenotype ranking captured growth,
muscle and perinatal features but did not explain the distinctive cancer signal
as well as BUB1B.

## Sensitivity and error analysis

Expanding phenotype retrieval from 500 to 2,000 genes introduced additional
pathogenic and rare candidates, including LZTR1 and GNRHR, but did not displace
the BUB1B pair.

An additional phenotype-weight ablation set all supplied HPO feature weights to
1.0 and changed missing population AF from a favorable score to a neutral score.
The BUB1B pair remained first with hypothesis score 0.9223, followed by the PEX5
homozygous deletion at 0.6992. This showed that the leading result was not
dependent on the initial clinical feature weights.

Important rejected or downgraded patterns included:

- ATP6V1E1 variants phased in cis;
- clustered SERPINA1 calls that generated false pair combinations;
- dense HLA-DRB1 haplotype calls;
- common variants with contradictory clinical labels;
- single heterozygous alleles in recessive-only models.

## Targeted raw-read validation

Lane-1 R1 reads were independently aligned to GRCh38 after ranking. At
`chr15:40209701`, three of four quality-filtered reads supported the BUB1B
alternate allele. At `chr15:40220612`, five of ten reads supported the alternate
allele. This provided independent raw-read confirmation that both submitted
BUB1B alleles were present in Lane 1.

The PEX5 deletion could not be reliably assessed using single-end R1 pileup.
Only two reference-supporting reads were visible in that lane, but ordinary
pileup is not sufficient to confirm or reject a 48-base deletion. PEX5 therefore
remained the strongest homozygous alternative based on the original multisource
VCF evidence.

Candidates below the leading hypotheses were not individually realigned. Their
ranking used VCF technical metrics, phenotype similarity, consequence,
population frequency, clinical evidence and inheritance compatibility.

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

## Innovation

The main methodological contribution is an uncertainty-preserving retrieval and
ranking architecture for an N-of-1 genome:

1. **Ranking rather than irreversible filtering.** The complete 5.0-million-row
   feature store was retained while separate candidate pools represented
   functional, homozygous, heterozygous-pair and mosaic hypotheses.
2. **Phenotype-cluster scoring.** HPO similarity was combined with coverage
   across cancer, renal, growth, muscle, perinatal and family-history groups.
   This reduced the risk that one common or nonspecific feature dominated gene
   retrieval.
3. **Hypothesis-aware ranking.** Single heterozygous, homozygous and
   compound-heterozygous models were scored separately rather than treating all
   variants as interchangeable rows.
4. **Pair-specific error controls.** The pipeline used phase fields, local
   distance, rarity of both alleles and inheritance compatibility to prevent
   cis variants, normalization fragments and dense haplotypes from becoming
   false compound-heterozygous predictions.
5. **Leaderboard-aligned evaluation.** The final representation directly
   modeled the challenge unit of prediction: a ranked single variant or
   two-variant pair with confidence thresholds.

## Scalability

The workflow was designed so that expensive operations act only on progressively
smaller datasets:

```text
5,012,204 complete variants
→ one reusable Parquet feature store
→ 6,390 Top-500 coding/splice candidates
→ 24,898 Top-2000 sensitivity candidates
→ tens of rare functional candidates
→ ranked hypotheses
→ targeted read validation
```

Scalability features include:

- streaming VCF projection with `bcftools`;
- lazy columnar feature engineering with Polars and Parquet;
- cached intermediate datasets to avoid rescanning the VCF;
- VEP annotation restricted to retrieved coding/splice regions;
- configurable HPO Top-K expansion;
- separate high-recall mosaic storage rather than expensive immediate analysis;
- targeted FASTQ/BAM validation only after candidate ranking;
- explicit row counts and sensitivity comparisons at each stage.

The same components can be reused for another proband by replacing the VCF,
reference-compatible annotations and HPO configuration. Thresholds and Top-K
scope remain configuration parameters rather than case-specific code changes.

## Work not claimed as completed

The current submission does not claim:

- a trained supervised or learning-to-rank model;
- confirmed parental segregation or BUB1B trans phase;
- completed CNV/SV event calling;
- complete paired-end FASTQ realignment.

These remain possible extensions. They were not presented as completed evidence
because the current ranking could be produced more efficiently and reproducibly
from the VCF, phenotype and targeted annotation data.

## Limitations

- Parental genomes were unavailable, so BUB1B trans configuration was not
  confirmed.
- R2 FASTQs were unavailable locally because of storage constraints.
- Only Lane-1 R1 reads were realigned; paired-end and cross-lane validation were
  not completed.
- Missing population frequency was treated as unknown, not automatically rare.
- HPO annotations incompletely represent family reproductive history and
  mechanism-level relationships.
- EPCR values are normalized hypothesis scores rather than formally calibrated
  probabilities.

## Reproducibility

The repository includes scripts for:

- VCF feature extraction;
- Parquet conversion;
- local HPO gene ranking;
- coding/splice region generation;
- targeted VEP annotation;
- structured annotation extraction;
- candidate and pair ranking;
- inheritance-aware hypothesis ranking;
- submission validation.

Public resources:

- Ensembl VEP release 116 GRCh38 cache
- Ensembl release 116 GRCh38 FASTA and GTF
- Human Phenotype Ontology `hp.obo`
- HPO `genes_to_phenotype.txt`

## AI assistance disclosure

GitHub Copilot CLI powered by GPT-5.6 Sol was used for pipeline design, code
generation, debugging, documentation and interpretation support. The account
plan/tier was not identified during the session. Raw VCF and FASTQ files were
processed locally and were not uploaded through the assistant. AI interaction
used aggregate statistics, selected candidate-level excerpts and public
challenge/reference information.
