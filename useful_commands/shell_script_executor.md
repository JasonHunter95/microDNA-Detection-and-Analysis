**microDNA-Detection-and-Analysis: Circle‑Map Pipeline Runner**

---

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Installation](#installation)
4. [Directory Structure](#directory-structure)
5. [Pipeline Scripts](#pipeline-scripts)
   - [Master Pipeline](#master-pipeline)
   - [Individual Scripts](#individual-scripts)
6. [Usage Examples](#usage-examples)
7. [Customizing & Adding New Scripts](#customizing--adding-new-scripts)
8. [Contributing](#contributing)
9. [License](#license)

---

## Overview

This repository implements a fully automated, reproducible pipeline for detecting, extracting, and annotating small extrachromosomal circular DNA (eccDNA) using [Circle‑Map](https://github.com/iprada/Circle-Map). The pipeline covers every step from reference download through final gene annotation, organized into modular scripts.

## Prerequisites

- **Operating System:** macOS or Linux
- **Software Dependencies:**
  - [BWA](http://bio-bwa.sourceforge.net/)
  - [Samtools](http://www.htslib.org/)
  - [Circle‑Map](https://github.com/iprada/Circle-Map)
  - [BEDTools](https://bedtools.readthedocs.io/)
  - Python ≥ 3.7 with required packages (`requirements.txt`)
  - `curl`, `gunzip`, `bash`
- Ensure all executables are in your `$PATH`.

## Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/YourUsername/microDNA-Detection-and-Analysis.git
   cd microDNA-Detection-and-Analysis
   ```
2. **Install Python dependencies** (if applicable)
   ```bash
   pip install -r requirements.txt
   ```
3. **Make scripts executable**
   ```bash
   chmod +x src/utils/shell_scripts/*.sh
   chmod +x src/pipeline/run_pipeline.sh
   ```

## Directory Structure

```
microDNA-Detection-and-Analysis/
├── data/                  # Reference, BAMs, BEDs, results
├── src/
│   ├── utils/
│   │   └── shell_scripts/  # Individual helper scripts
│   └── pipeline/           # Master pipeline script
├── requirements.txt       # Python dependencies
└── README.md
```

## Pipeline Scripts

### Master Pipeline

The master runner orchestrates all steps across chromosomes:

```bash
bash src/pipeline/run_pipeline.sh \
  -s <SAMPLE_ID>  \
  -t <THREADS>     \
  -o <OUTPUT_DIR>
```

**Options:**

- `-s`: Sample ID (e.g., SRR413984)
- `-t`: Number of threads (default: 8)
- `-o`: Base output directory (default: `data`)

### Individual Scripts

Each shell script can also be run independently.

#### 1. Download Reference Genome

```bash
bash src/utils/shell_scripts/get_ref_genome.sh
```

#### 2. Extract Specific Chromosome

```bash
bash src/utils/shell_scripts/extract_chr.sh -c <CHR>
# e.g., -c 1 extracts chr1
```

#### 3. Index All Chromosomes

```bash
bash src/utils/shell_scripts/index_all_chr.sh
```

#### 4. Align Reads & Sort BAM by Chromosome

```bash
bash src/utils/shell_scripts/bwa_align_and_sort_by_chr.sh \
  <CHR> <SAMPLE_ID>
```

#### 5. Query-Sort BAM

```bash
bash src/utils/shell_scripts/querysort_bam_by_chr.sh \
  <CHR> <SAMPLE_ID>
```

#### 6. Sort & Index Candidate BAMs

```bash
bash src/utils/shell_scripts/sort_and_index_candidates_by_chr.sh \
  <CHR> <SAMPLE_ID>
```

#### 7. Clean BED Output

```bash
python3 src/utils/clean_bed.py \
  -i data/beds/chr${CHR}/${SAMPLE_ID}_chr${CHR}.eccdna.bed \
  -o data/beds/chr${CHR}/${SAMPLE_ID}_chr${CHR}.eccdna.cleaned.bed
```

#### 8. Gene Annotation

```bash
gtf2bed < data/gtfs/gencode.annotation.gtf > annotation.bed
bedtools intersect -a cleaned.bed -b annotation.bed -wa -wb > with_genes.bed
```

## Usage Examples

1. **Run full pipeline for sample SRR413984 on 4 threads**
   ```bash
   bash src/pipeline/run_pipeline.sh -s SRR413984 -t 4 -o data
   ```
2. **Extract only chromosome X for quick test**
   ```bash
   bash src/utils/shell_scripts/extract_chr.sh -c X
   ```

## Customizing & Adding New Scripts

1. **Add your script** to `src/utils/shell_scripts/` or `src/pipeline/`.
2. **Grant execute permissions**:
   ```bash
   chmod +x src/utils/shell_scripts/your_script.sh
   ```
3. **Document** the new script under [Individual Scripts](#individual-scripts) with usage examples.
4. **Optionally**, integrate into `run_pipeline.sh` by adding a step in the `for chr` loop.

## Contributing

- Fork this repository
- Create a feature branch (`git checkout -b feature/my-script`)
- Commit your changes (`git commit -m "Add my new script"`)
- Push to your fork (`git push origin feature/my-script`)
- Open a Pull Request

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.