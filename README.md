# microDNA Detection and Analysis

[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-32%20passed-success.svg)]()

A reproducible bioinformatics pipeline for detecting, extracting, and annotating **microDNA** (small extrachromosomal circular DNA, 100–400 bp) from high-throughput sequencing data.

---

## Overview

Extrachromosomal circular DNA (eccDNA) are DNA molecules found outside of chromosomes that have been implicated in gene amplification, aging, and cancer. This pipeline leverages [Circle-Map](https://github.com/iprada/Circle-Map) to identify eccDNA candidates from paired-end sequencing reads, then cleans, filters, and annotates them with genomic features from GENCODE.

### Key Features

- **Automated Detection**: Uses Circle-Map's statistical approach to identify circular DNA candidates.
- **Size Filtering**: Filters candidates to the microDNA size range (100–400 bp).
- **Gene Annotation**: Annotates eccDNA with nearest genes using GENCODE.
- **Reproducible**: Conda environment ensures consistent dependencies.
- **Modular CLI Scripts**: Each pipeline step is a separate, reusable command-line tool.

---

## Project Structure

```
microDNA-Detection-and-Analysis/
├── data/                       # Data directory (created during setup)
│   ├── fastqfiles/             # Raw FASTQ reads
│   ├── human_genome/           # Reference genome files
│   ├── bams/                   # Alignment files
│   ├── beds/                   # BED files (eccDNA coordinates)
│   └── results/                # Final annotated outputs
├── src/utils/
│   ├── py_scripts/             # Python CLI scripts
│   │   ├── microDNA.py         # Soft-clip based microDNA detection
│   │   ├── clean_bed.py        # BED formatting and filtering
│   │   ├── annotate_eccdna_closest.py    # Gene annotation (closest)
│   │   ├── annotate_eccdna_intersected.py # Gene annotation (overlap)
│   │   ├── gtf_parser.py       # Shared GTF/GENCODE parsing
│   │   ├── plotting.py         # Shared visualization utilities
│   │   └── ...                 # Additional utilities
│   └── shell_scripts/          # Shell utility scripts
├── tests/                      # Unit tests (32 tests)
├── environment.yml             # Conda environment definition
├── requirements.txt            # Python dependencies
└── README.md
```

---

## Quick Start

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)
- ~50 GB disk space (for reference genome and example data)

### Installation

```bash
# Clone the repository
git clone https://github.com/JasonHunter95/microDNA-Detection-and-Analysis.git
cd microDNA-Detection-and-Analysis

# Create and activate the conda environment
conda env create -f environment.yml
conda activate circlemap-env

# Install Python dependencies and the package in editable mode
pip install -e .
```

---

## Pipeline Workflow

```mermaid
flowchart LR
    A[FASTQ Reads] --> B[BWA-MEM2 Alignment]
    B --> C[Sorted BAM]
    C --> D[Query-sorted BAM]
    D --> E[Circle-Map ReadExtractor]
    E --> F[Candidate BAM]
    F --> G[Circle-Map Realign]
    G --> H[Raw eccDNA BED]
    H --> I[microdna.clean_bed]
    I --> J[bedtools closest]
    J --> K[microdna.annotate_eccdna_closest]
    K --> L[Annotated TSV]
```

---

## CLI Scripts Reference

All scripts are available as `microdna` modules.

| Script | Purpose |
|--------|---------|
| `clean_bed` | Format BED files to BED6 with unique eccDNA IDs |
| `annotate_eccdna_closest` | Annotate eccDNA with closest genes |
| `annotate_eccdna_intersected` | Annotate eccDNA with overlapping genes |
| `parse_intersected` | Parse bedtools intersect output |
| `unique_eccdna_table` | Aggregate eccDNA by unique ID |
| `cleaner_tsv` | Filter and validate eccDNA annotations |
| `eccdna_length_distribution` | Plot eccDNA length distribution |
| `eccdna_score_distribution` | Plot eccDNA score distribution |
| `get_soft_clips` | Count soft-clipped reads in BAM |
| `print_soft_clipped_reads` | Display soft-clipped read details |
| `sw` | Smith-Waterman local alignment |
| `get_seq` | Analyze eccDNA junction sequences |

---

## Detailed Usage

### 1. Data Acquisition

```bash
# Download sequencing reads
fasterq-dump SRR413984 --split-files -O data/fastqfiles/

# Download reference genome
bash src/utils/shell_scripts/get_ref_genome.sh

# Prepare GENCODE annotations
bash src/utils/shell_scripts/prepare_gencode_v19_annotation.sh
```

### 2. Alignment

```bash
# Index reference
bwa-mem2 index data/human_genome/chr1/chr1.fna
samtools faidx data/human_genome/chr1/chr1.fna

# Align and sort
bwa-mem2 mem -t 8 data/human_genome/chr1/chr1.fna \
    data/fastqfiles/SRR413984_1.fastq data/fastqfiles/SRR413984_2.fastq | \
    samtools sort -@ 8 -m 1G -o data/bams/chr1/SRR413984_chr1.sorted.bam

# Index BAM
samtools index data/bams/chr1/SRR413984_chr1.sorted.bam

# Query-sort for Circle-Map
samtools sort -n -o data/bams/chr1/SRR413984_chr1.querysorted.bam \
    data/bams/chr1/SRR413984_chr1.sorted.bam
```

### 3. eccDNA Detection

```bash
# Extract candidate circular reads
Circle-Map ReadExtractor \
    -i data/bams/chr1/SRR413984_chr1.querysorted.bam \
    -o data/bams/chr1/SRR413984_chr1.candidates.bam

# Sort and index candidates
samtools sort -o data/bams/chr1/SRR413984_chr1.candidates.sorted.bam \
    data/bams/chr1/SRR413984_chr1.candidates.bam
samtools index data/bams/chr1/SRR413984_chr1.candidates.sorted.bam

# Realign to detect eccDNA
Circle-Map Realign \
    -i data/bams/chr1/SRR413984_chr1.candidates.sorted.bam \
    -qbam data/bams/chr1/SRR413984_chr1.querysorted.bam \
    -sbam data/bams/chr1/SRR413984_chr1.sorted.bam \
    -fasta data/human_genome/chr1/chr1.fna \
    -o data/beds/chr1/SRR413984_chr1.eccdna.bed \
    --split 2 --threads 8
```

### 4. Post-processing

```bash
# Format to BED6 with unique IDs
python -m microdna.clean_bed \
    -i data/beds/chr1/SRR413984_chr1.eccdna.bed \
    -o data/beds/chr1/SRR413984_chr1.eccdna.cleaned.bed
```

### 5. Annotation

```bash
# Convert to UCSC chromosome naming
sed 's/^NC_000001.10/chr1/' data/beds/chr1/SRR413984_chr1.eccdna.cleaned.bed \
    > data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.bed

# Sort eccDNA BED
bedtools sort -i data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.bed \
    > data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.sorted.bed

# Find closest genes
bedtools closest \
    -a data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.sorted.bed \
    -b data/beds/whole_genomes/human_genome/gencode.v19.annotation.genes.sorted.bed \
    -d > data/beds/chr1/SRR413984_chr1.eccdna.closest_genes.bed

# Generate annotated output
python -m microdna.annotate_eccdna_closest \
    -i data/beds/chr1/SRR413984_chr1.eccdna.closest_genes.bed \
    -o data/results/chr1/SRR413984_chr1.eccdna.annotated.tsv
```

---

## Output Format

The final annotated TSV contains:

| Column | Description |
|--------|-------------|
| `ecc_id` | Unique eccDNA identifier (e.g., `ecc_00001`) |
| `ecc_chrom` | Chromosome |
| `ecc_start` | Start position (0-based) |
| `ecc_end` | End position |
| `gene_id` | ENSEMBL gene ID |
| `gene_name` | Gene symbol |
| `gene_type` | Gene biotype (protein_coding, lncRNA, etc.) |
| `distance` | Distance to nearest gene (0 if overlapping) |

---

## Example Results

### eccDNA Visualization in IGV

Circular DNA candidates detected on chromosome 1, visualized in the Integrative Genomics Viewer:

![eccDNA detected on chr1](figures/eccDNA_on_chr1.png)

### Distribution Analysis

The pipeline generates distribution plots to characterize detected eccDNA:

| Length Distribution | Score Distribution |
|:---:|:---:|
| ![eccDNA length distribution](figures/eccDNA_length_distribution.png) | ![eccDNA score distribution](figures/eccDNA_score_distribution.png) |

---

## Testing

Run the test suite:

```bash
conda activate circlemap-env
pytest tests/ -v -n auto
```

**Current status**: 32 tests passing

---

## Development

### Linting

```bash
ruff check src/
```

### Code Style

This project uses:
- **ruff** for linting
- **Type hints** for function signatures
- **Docstrings** for all public functions
- **argparse** for CLI interfaces

---

## Windows Setup (WSL)

This pipeline requires a Unix environment. On Windows, use **Windows Subsystem for Linux (WSL2)**:

1. **Install WSL2** (PowerShell as Administrator):
   ```powershell
   wsl --install
   ```

2. **After restart**, open Ubuntu and install Miniconda:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
   bash miniconda.sh
   source ~/.bashrc
   ```

3. **Clone and set up** (from WSL terminal):
   ```bash
   cd ~
   git clone https://github.com/JasonHunter95/microDNA-Detection-and-Analysis.git
   cd microDNA-Detection-and-Analysis
   conda env create -f environment.yml
   conda activate circlemap-env
   ```

> **Tip**: For better performance, work within the WSL filesystem (`~/`) rather than `/mnt/c/`.

---

## References

- **Circle-Map**: Prada-Luengo, I., Krogh, A., Maretty, L., & Regenberg, B. (2019). Sensitive detection of circular DNAs at single-nucleotide resolution using guided realignment of partially aligned reads. *BMC Bioinformatics*, 20(1), 1-9. [GitHub](https://github.com/iprada/Circle-Map)
- **GENCODE**: [gencodegenes.org](https://www.gencodegenes.org/)

---
