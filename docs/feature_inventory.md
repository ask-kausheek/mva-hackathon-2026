# B. Feature inventory

Features are grouped by the ranking signal they provide.

| Group | Features | Ranking purpose | Cost |
|---|---|---|---|
| Technical | QUAL, FILTER, DP, GQ, AD, VAF, MQ, QD, FS, SOR, rank-sum metrics | Remove obvious artifacts; preserve low-VAF evidence separately | Low |
| Variant | SNP/indel, allele lengths, multiallelic, rsID, GT, PGT/PID phase set | Route variants to appropriate models and pair candidates | Low |
| Population | gnomAD AF/AC/AN, homozygotes, ancestry max AF | Strong rare-variant retrieval | Medium, local join |
| Functional | gene, transcript, consequence, coding/splice class, protein change | Reduce millions to functionally plausible candidates | Medium |
| Pathogenicity | ClinVar, CADD/REVEL/AlphaMissense/SpliceAI when locally available | Prioritize previously supported or damaging variants | Medium |
| Gene-disease | inheritance mode, disease validity, cancer/MVA association | Constrain valid compound-het pairs | Medium |
| Phenotype | HPO gene score, matched terms, information content | Case-specific ranking rather than generic pathogenicity | Medium |
| Pair | same gene, two heterozygous variants, phase evidence, consequence combination | Test the organizer-signaled compound-het hypothesis | Low after retrieval |
| MVA/CNV | VAF deviation, chromosome depth, B-allele imbalance, breakpoint/SV evidence | Separate mosaic/chromosome branch | High; defer |

## Feature policy

1. Technical quality is a reliability feature, not the causal target.
2. Population rarity is evaluated from an external population database, never call-set `INFO/AF`.
3. Hard filters are limited to clear invalid records. Borderline records remain in alternate pools.
4. Gene-level phenotype scores are joined once and reused across all variants in that gene.
5. Pair features are computed only after single-variant retrieval has reduced the search space.
6. PASS, dbSNP membership, and high QUAL are never treated as stand-alone pathogenicity evidence.

## Required local annotation minimum

The first useful annotation pass needs only:

- normalized GRCh38 key: chromosome, position, reference, alternate
- gene
- consequence
- population maximum allele frequency
- ClinVar significance

VEP is optional. Any local annotator or prebuilt GRCh38 table is acceptable if it supplies these features reproducibly.

For the installed VEP cache, the targeted pass requests:

- `SYMBOL`, consequence, biotype and variant class
- canonical/MANE transcript indicators
- coding and protein HGVS
- 1000 Genomes and gnomAD allele frequencies when present in the cache
- maximum colocated population frequency
- clinical significance when present

Missing cache fields remain null. Missing population AF receives a neutral score,
not a rarity bonus. Missing clinical evidence is not interpreted as benign.

## Implemented score policy

- Phenotype relevance, consequence, rarity, clinical evidence and technical
  evidence are normalized before weighted combination.
- Clinical-prior HPO weights were compared with an all-equal HPO-weight
  sensitivity run.
- Conflicting pathogenic/benign labels do not receive a pathogenic override.
- Single-variant hypotheses are gated by inheritance compatibility, clinical
  classification and phenotype rank.
- Final EPCR values are a monotonic min-max transformation of the selected
  Top-10 hypothesis scores, not calibrated probabilities.
