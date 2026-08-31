# Inheritance Hypothesis Testing

## What we actually know

The compound-heterozygous statement came from the public Track 1
`evaluation.py` docstring:

> simplified for a single proband ("N-of-1") with a clinically validated
> compound-heterozygous answer key

The scorer also contains special partial-credit logic for recovering one of two
true variants. This is direct organizer-authored metadata about the production
answer structure.

However:

- the actual production ground-truth coordinates are private;
- `groundtruth.py` exposes only fake local fallback coordinates;
- we have not independently demonstrated two causal variants from the VCF;
- the 11 positive-control candidates were selected using a heterozygous filter, so that
  subset cannot fairly compare all inheritance models.

Therefore, compound heterozygosity is a **high-priority hypothesis and submission
format clue**, not a conclusion that replaces data analysis.

## Competing hypotheses

| ID | Hypothesis | Candidate evidence | Evidence against |
|---|---|---|---|
| H1 | Compound heterozygous | Two rare damaging heterozygous variants in one relevant gene, preferably in trans | Common/benign variants, same haplotype, only one credible allele |
| H2 | Homozygous recessive | Rare damaging `1/1` variant in a relevant gene with credible depth | Common homozygous variant or weak technical support |
| H3 | Dominant/de novo-like single variant | One rare high-impact heterozygous variant with strong phenotype fit | No parental data, common population frequency, weak consequence |
| H4 | Mosaic SNV/indel | Reproducible low VAF with strong read support and relevant consequence | Ordinary allele imbalance, mapping artifacts, dbSNP/common variant |
| H5 | CNV/SV/aneuploidy | Depth, allele-balance or breakpoint signal matching the phenotype | No reproducible event-level signal |

Without parental genomes, “de novo” and “in trans” cannot be proven from the
proband VCF alone. Read phasing can support trans/cis locally; parental or
long-range evidence is stronger.

## Fair testing design

Do not compare hypotheses using the current 11-row table alone. That table was
created by:

```text
three selected genes + heterozygous genotype + technical thresholds
```

It is appropriate for testing H1 inside those genes, but it excludes H2 by
construction and does not test H3-H5 genome-wide.

Use separate candidate pools:

```text
Pool A: rare damaging heterozygous variants
Pool B: rare damaging homozygous variants
Pool C: low-VAF candidates
Pool D: CNV/SV/chromosome events
Pool E: known pathogenic variants regardless of inheritance model
```

Rank all pools using comparable feature groups:

- population rarity;
- functional consequence;
- clinical evidence;
- phenotype relevance;
- technical evidence;
- inheritance compatibility.

## Phenotype-driven retrieval

The eight supplied HPO terms are converted into a local gene ranking using the
public HPO ontology and gene annotations. The baseline:

1. propagates gene annotations through HPO ancestor relationships;
2. calculates information-content semantic similarity;
3. weights the supplied clinical features;
4. rewards coverage across cancer, renal, growth, muscle, perinatal and family
   history groups.

Recurrent spontaneous abortion is retained as a parental/family-history feature.
The group-coverage term prevents a gene matching only one feature from dominating
the complete clinical constellation.

## Positive-control result

The positive-control analysis asked:

1. Do the literature-derived genes `BUB1B`, `CEP57`, or `TRIP13` contain two
   rare functional heterozygous variants?
2. Are those variants technically credible?
3. Is phase evidence compatible with trans configuration?

It identified the two BUB1B variants later recovered as the leading pair in both
the independent Top-500 and Top-2000 phenotype-driven analyses. The hypothesis
was therefore retained based on converging evidence rather than the seed list
alone.

## Evidence standard before submission

A pair should not rank first merely because the scorer can accept two variants.
It should have:

- two independently credible calls;
- low population frequency;
- relevant functional consequences;
- same-gene and inheritance compatibility;
- phenotype-relevant gene evidence;
- no evidence that both variants are in cis;
- targeted read validation for final candidates.

The final Top-10 can contain a mixture of pair and single-variant rows if the
evidence remains uncertain.

## Final sensitivity result

The phenotype analysis was repeated with all HPO feature weights equal to 1.0
and missing population frequency treated as neutral. BUB1B remained the leading
pair, followed by the PEX5 homozygous deletion. This reduced concern that the
result was created by the initial clinical weighting scheme.

Lane-1 R1 alignment independently observed both BUB1B alternate alleles. Phase
remains unknown, and the R1 validation does not replace parental segregation.
