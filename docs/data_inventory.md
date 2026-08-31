# A. Data inventory

## Available case data

| Asset | Role | Current decision |
|---|---|---|
| `WGS_EX2312012_HGWCNDSX7.vcf.gz` | Primary candidate source | Process first |
| VCF tabix index | Random region access | Keep beside VCF |
| Four local R1 FASTQ lanes | Read-level confirmation | Lane 1 used to confirm both BUB1B alleles |
| Eight HPO terms | Phenotype retrieval and gene ranking | Use immediately |

## VCF profile

| Statistic | Value |
|---|---:|
| Samples | 1 |
| Records | 5,012,204 |
| SNPs | 4,086,641 |
| Indels | 930,181 |
| Multiallelic sites | 66,860 |
| PASS | 4,740,790 |
| `0/1` | 3,047,077 |
| `1/1` | 1,898,267 |
| `1/2` | 66,860 |

The VCF contains technical calling fields. No provided pre-annotated variant table, population resource, or clinical annotation resource has been identified. `INFO/AF` is sample/cohort call-set frequency here, not a population-frequency feature; it must not be used as rarity evidence.

## Observed record patterns

- VCF chromosomes are represented as `1`, while the submission template demonstrates `chr1`. The feature table therefore stores both the raw chromosome and a canonical `submission_chrom`.
- PASS is not a sufficient quality filter. An observed PASS homozygous call has depth 2 and GQ 6, so DP, GQ, AD/VAF, and mapping metrics remain ranking features.
- dbSNP identifiers are common. An `rs` identifier means the variant is catalogued, not benign or common; it is retained as a weak annotation feature only.
- Nearby heterozygous variants can share `PID` and phased genotype `PGT`. These local phase blocks are useful for identifying variants on the same haplotype and avoiding biologically invalid compound-heterozygous pairs.
- Indels are present among high-quality PASS calls and must be normalized before matching external annotation databases or constructing submission coordinates.

## Public challenge resources

- Submission CSV template
- Methods-description spreadsheet
- Public scoring implementation
- Official rules, FAQ, and upload UI

The public scorer's docstring states that its production answer key is clinically
validated and compound heterozygous. This is organizer metadata, not independent
verification from our VCF. A full pair match requires both variants in one row;
single-variant rows are also supported. Coordinates must be GRCh38.

## Submission schema

`proband_id, chrom_1, pos_1, ref_1, alt_1, chrom_2, pos_2, ref_2, alt_2, epcr, finding_type, notes`

Constraints:

- `proband_id` is `PROBAND01`
- At most 10 rows
- `epcr` is in `(0, 1]`; rows are sorted descending by it
- Six submissions per participant
- Report is PDF or Markdown
- A GitHub URL is required; the repository may remain private during the
  hackathon and must be public when judging requires it

## Data-governance constraint

Raw and derived individual-level data cannot be redistributed and must be deleted within the rule-defined retention period. Candidate tables therefore stay under ignored `data/` and `results/` paths.
