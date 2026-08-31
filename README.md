# SageBio MVA Hackathon 2026

Data-engineering and ranking pipeline for Track 1. Raw genomic and clinical data must never be committed.

## Submission status

Submission 1 is prepared:

```text
submission/ask-kausheek_hpo-vep-rankfusion.csv
submission/ask-kausheek_track1_report.md
```

## Confirmed inputs

- `<WGS_PROBAND>.vcf.gz` and its tabix index
- One sample, GRCh38, 5,012,204 records
- FASTQ files reserved for targeted validation
- Eight supplied HPO terms in `configs/proband.yaml`
- Official submission limit: 10 ranked rows and 6 submissions

## First run

From WSL, with the VCF in the current directory:

```bash
export WGS_PROBAND_VCF=/path/to/proband.vcf.gz

./scripts/01_extract_vcf_features.sh "$WGS_PROBAND_VCF" data/interim
python -m src.features.build_vcf_features \
  --input data/interim/vcf_features.tsv.gz \
  --output-dir data/processed
```

This performs one full VCF scan, stores reusable columnar features, and creates a separate mosaic-review pool. It does not annotate variants or process FASTQs.

## Targeted VEP pass

VEP should first run on prioritized gene regions rather than all 5 million records:

```bash
./scripts/02_targeted_vep.sh \
  "$WGS_PROBAND_VCF" \
  /path/to/GRCh38.fa \
  configs/known_mva_seed_genes.txt \
  /path/to/Homo_sapiens.GRCh38.gtf.gz \
  /path/to/vep_data \
  data/annotation/tier1 \
  8
```

The script extracts gene regions before VEP, keeps credible heterozygous calls, normalizes alleles, and requests only ranking-relevant annotations.

For the Docker installation, use the split workflow:

```bash
./scripts/02_prepare_targeted_candidates.sh \
  "$WGS_PROBAND_VCF" GRCh38.fa \
  configs/known_mva_seed_genes.txt Homo_sapiens.GRCh38.gtf.gz \
  data/annotation/tier1

./scripts/03_run_vep_docker.sh \
  GRCh38.fa vep_data data/annotation/tier1 8
```

## Design documents

- [System architecture](ARCHITECTURE.md)
- [Inheritance hypothesis testing](HYPOTHESIS_TESTING.md)
- [Observed results and ranking evidence](OBSERVATIONS.md)
- [Chronological session summary](SESSION.md)
- [Data inventory](docs/data_inventory.md)
- [Feature inventory](docs/feature_inventory.md)
- [Pipeline design](docs/pipeline_design.md)
- [First baseline](docs/first_baseline.md)
- [Experiment plan](docs/experiment_plan.md)

## Privacy

The gated data, sample-level derivatives, candidate tables, and reports containing private clinical/genomic details remain local and must be deleted under the hackathon rules. Only reusable code, aggregate statistics, and approved submission artifacts belong in GitHub.
