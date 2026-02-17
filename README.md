# STREAM 🌊

**S**treamlined **T**ranscript **E**xpression & **R**NA-seq **M**apping

Nextflow DSL2 pipeline for RNA-seq quality control and transcript quantification.
Ultra-minimalist — designed for solo bioinformaticians. Inspired by [nf-core/rnaseq](https://nf-co.re/rnaseq).

---

## Pipeline Overview

```
  FASTQ / SRA / GEO
         │
         ▼
  ┌─ 1. FastQC (raw) ───────────────────────────────────────┐
  │                                                          │
  ├─ 2. fastp (trimming + QC) ──────────────────────────────┤
  │         │                                                │
  │         ├─ 3. FastQC (clean) ───────────────────────────┤
  │         ├─ 4. FastQ Screen (opt) ──────────────────────-┤
  │         ├─ 5. Sequence stats ──────────────────────────-┤
  │         ├─ 6. Kraken2 (opt) ───────────────────────────-┤
  │         └─ 7→8. Salmon index + quant (opt) ────────────-┤
  │                                                          │
  └──────────────── 9. MultiQC (aggregation) ───────────────┘
```

## Quick Start

```bash
# From a FASTQ directory (auto-detects PE/SE)
nextflow run BIPN/STREAM --fastq_dir /path/to/fastqs --outdir results -resume

# From a samplesheet CSV
nextflow run BIPN/STREAM --input samplesheet.csv --outdir results -resume

# From SRA accessions
nextflow run BIPN/STREAM --sra_ids "SRR1234567,SRR1234568" --outdir results -resume

# From a GEO dataset (auto-resolves GSE → SRR)
nextflow run BIPN/STREAM --sra_ids GSE123456 --outdir results -resume

# Full pipeline with all QC options
nextflow run BIPN/STREAM \
    --fastq_dir /path/to/fastqs \
    --run_salmon \
    --run_fastq_screen --fastq_screen_conf /path/to/fastq_screen.conf \
    --run_kraken2 --kraken2_db /path/to/kraken2_db \
    --outdir results \
    -resume
```

### Samplesheet format (CSV)

```csv
sample,fastq_1,fastq_2
sampleA,/path/to/sampleA_R1_001.fastq.gz,/path/to/sampleA_R2_001.fastq.gz
sampleB,/path/to/sampleB_R1_001.fastq.gz,
```

> `fastq_2` vide = single-end. Les fichiers multi-lanes avec le même `sample` sont traités séparément. Pour merger les lanes, pré-concaténer ou dupliquer les lignes dans le samplesheet.

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `null` | Samplesheet CSV (`sample,fastq_1,fastq_2`) |
| `--fastq_dir` | `null` | Dossier de FASTQs (`*_R{1,2}_001.fastq.gz`) |
| `--sra_ids` | `null` | Accessions SRA/GEO (CSV ou fichier, 1/ligne) |
| `--outdir` | `results` | Dossier de sortie |
| `--run_salmon` | `true` | Activer quantification Salmon |
| `--salmon_index` | `null` | Index Salmon pré-construit |
| `--transcriptome_fasta` | `null` | FASTA transcriptome (skip download) |
| `--genome` | `GRCh38` | Assemblage génomique (Ensembl) |
| `--ensembl_release` | `115` | Version Ensembl |
| `--run_fastq_screen` | `false` | Activer FastQ Screen |
| `--fastq_screen_conf` | `null` | Fichier config FastQ Screen |
| `--run_kraken2` | `false` | Activer Kraken2 |
| `--kraken2_db` | `null` | Base Kraken2 |
| `--fastp_qualified_quality` | `20` | Phred score minimum (fastp) |
| `--fastp_length_required` | `20` | Longueur minimum après trim |
| `--skip_fastqc` | `false` | Désactiver FastQC |
| `--save_trimmed` | `false` | Publier les FASTQs trimmés |
| `--subset_size` | `0` | FastQ Screen subset (0 = all) |
| `--max_cpus` | auto | Nombre max de CPUs |

---

## Output Structure

```
results/
├── 00_sra_fastq/        # FASTQs téléchargés (si SRA)
├── 01_fastqc_raw/       # QC reads bruts
├── 02_fastp/            # Rapports trimming + FASTQs (si --save_trimmed)
├── 03_fastqc_clean/     # QC post-trimming
├── 04_fastq_screen/     # Contamination screening (opt)
├── 05_statistics/       # Stats séquences (seqtk)
├── 06_kraken2/          # Classification taxonomique (opt)
├── 07_salmon/           # Quantification transcripts (opt)
├── 08_multiqc/          # Rapport agrégé interactif
├── reference/           # Transcriptome + index Salmon (cache)
└── pipeline_info/       # Nextflow timeline, trace, DAG, report
```

---

## Requirements

**Core** (toujours requis) :
`fastqc` `fastp` `multiqc` `seqtk`

**Optionnel** :
`salmon` (quantification) · `fastq_screen` `bowtie2` (contamination) · `kraken2` (taxonomie) · `sra-tools` `pigz` (SRA download)

**Nextflow** ≥ 23.04

---

## Resume & Cache

Le pipeline exploite nativement le cache Nextflow (`-resume`). Les étapes déjà complétées sont automatiquement sautées. Les références (transcriptome, index Salmon) sont persistées via `storeDir` et réutilisées entre exécutions.

```bash
# Relancer après un crash — reprend exactement là où ça s'est arrêté
nextflow run main.nf --fastq_dir fastqs --outdir results -resume
```

---

## License

MIT
