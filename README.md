# microDNA-Detection-and-Analysis

**A reproducible pipeline for detecting, extracting, and annotating small extrachromosomal circular DNA (eccDNA) from high‑throughput sequencing data.**

---

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)  
3. [Installation](#installation)  
4. [Data Acquisition](#data-acquisition)  
5. [Alignment](#alignment)  
6. [eccDNA Detection](#eccdna-detection)  
7. [Post‑processing](#post-processing)  
8. [Annotation](#annotation)
9. [Results](#results)

---

## Overview

This pipeline uses [Circle‑Map](https://github.com/iprada/Circle-Map) to detect eccDNA (100–400 bp) from paired‑end sequencing reads, cleans and re‑aligns candidates, and annotates genomic features. It is modular, transparent, and version‑controlled for full reproducibility.
It is designed to be run on a Linux system with a conda package manager.

---

## Prerequisites

- **Conda** (miniconda or anaconda)
- Python (3.6 or later)
- [BWA](http://bio-bwa.sourceforge.net/)  
- [Samtools](http://www.htslib.org/)  
- [Bedtools](https://bedtools.readthedocs.io/)  
- [gtf2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gtf2bed.html)

## Installation

1. Clone the Repository

    ```bash
    git clone https://github.com/JasonHunter95/microDNA-Detection-and-Analysis.git
    cd microDNA-Detection-and-Analysis
    ```

2. Create and Activate a Conda Environment

    ```bash
    conda env create -f environment.yml
    conda activate circlemap-env
    ```

---

## Data Acquisition

To run the pipeline, you need to have paired-end sequencing data in FASTQ format. You can either download your own data or use the example dataset provided in this repository.

Raw sequencing data can be fetched from the SRA database. The example dataset is from the NCBI SRA database, specifically the SRR413969 accession number.

```bash
fasterq-dump SRR413969 --split-files -O data/
```

This will download the paired-end reads and save them in the `data/` directory. The files will be named `SRR413969_1.fastq` and `SRR413969_2.fastq`.

The example dataset is a human sample, and the reference genome used in this example is the human genome (GRCh37). The reference genome file can be downloaded to the `data/` directory as `GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna` via executing the following shell script:

Please note that this may take some time to download, as the reference genome is quite large.

```bash
bash src/utils/shell_scripts/get_ref_genome_1.sh
```

## Alignment

The first step in the pipeline is to align the paired-end reads to the reference genome. The reference genome used in this example is the human genome (GRCh37).
This may take some time to index, please be patient.

```bash
bwa index data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna
```

Align the reads to the reference genome using BWA, and then sort the output BAM file using Samtools. The file is sent to the data/ folder. The `-t` option specifies the threadcount to use for BWA. I only used 4 threads for this because it was done locally on my mac. This can be increased if you are running on a more powerful machine.

```bash
bwa mem -t 4 data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna \
    data/SRR413969_1.fastq data/SRR413969_2.fastq | \
    samtools sort -o data/SRR413969.sorted.bam
```

Now, index the sorted BAM file using Samtools:

```bash
samtools index data/SRR413969.sorted.bam
```

Next, in order to prepare the BAM file for Circle-Map, we need to query sort the BAM file. This is done using the `samtools sort -n` command. The `-n` option sorts the reads by name, which is required for Circle-Map.

```bash
samtools sort -n -o data/SRR413969.querysorted.bam data/SRR413969.sorted.bam
```

## eccDNA Detection

The next step is to detect eccDNA using the Circle-Map tool. The Circle-Map pipeline consists of two main steps: ReadExtractor and Realign.
The ReadExtractor step extracts the reads that are most statistically likely to be circular, and the Realign step re-aligns the extracted reads to the reference genome.

```bash
Circle-Map ReadExtractor -i data/SRR413969.querysorted.bam -o data/SRR413969.candidates.bam
```

Then sort and index the candidates.bam file using Samtools:

```bash
samtools sort -o data/SRR413969.candidates.sorted.bam data/SRR413969.candidates.bam
samtools index data/SRR413969.candidates.sorted.bam
```

Now we can run the Realign step. This step takes the queryname sorted BAM file and the original sorted BAM file as input, and outputs a BED file with the eccDNA coordinates.

```bash
Circle-Map Realign \         
  -i ./data/SRR413969.candidates.sorted.bam \
  -qbam ./data/SRR413969.querysorted.bam \
  -sbam ./data/SRR413969.sorted.bam \
  -fasta ./data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna \
  -o ./data/SRR413969.eccdna.bed \
  --split 1 \
  --verbose 3
```

## Post-processing

Run the following script to clean the output files:

```bash
python3 src/utils/clean_bed.py \
  -i temp_files_<replace_with_your_numbers>/peaks.bed \
  -o data/SRR413969.eccdna.cleaned.bed
```

## Annotation

  1. Annotate the eccDNA using the Gencode annotation file.

      ```bash
        curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
        gunzip gencode.v19.annotation.gtf.gz
      ```

  2. Convert the GTF file to BED format using the `gtf2bed` command. This creates a BED file with the Gencode annotation.

      ```bash
        gtf2bed < gencode.v19.annotation.gtf > gencode.v19.annotation.bed
      ```

  3. Convert the eccDNA BED file to UCSC format. This is necessary for the `bedtools intersect` command to work. The Gencode annotation uses `chr1`, `chr2`, etc. as chromosome names, while the CircleMap output uses `NC_000001.10`, `NC_000002.11`, etc. as chromosome names.

      ```bash
        sed 's/^NC_000001.10/chr1/' data/SRR413969.eccdna.cleaned.bed > data/SRR413969.eccdna.cleaned.ucsc.bed
      ```

  4. Intersect the eccDNA BED file with the Gencode annotation BED file using `bedtools intersect`. This will create a new BED file with the eccDNA coordinates and their corresponding gene annotations.

      ```bash
        bedtools intersect -a data/SRR413969.eccdna.cleaned.ucsc.bed \
                          -b data/gencode.v19.annotation.bed \
                          -wa -wb > data/SRR413969.eccdna.with_genes.bed
      ```

  5. Generate a clean output of the eccDNA with gene annotations using the `annotate_eccdna_genes.py` script. This will create a new TSV file with the eccDNA coordinates and their corresponding gene annotations.

      ```bash
        python3 src/utils/annotate_eccdna_genes.py
      ```

## Results

The final output file `data/SRR413969.eccdna.annotated.tsv` contains the eccDNA coordinates and their corresponding gene annotations. The file is tab-separated and contains the following columns:

- `ecc_id`: eccDNA ID
- `chrom`: Chromosome name
- `start`: Start position of the eccDNA
- `end`: End position of the eccDNA
- `gene_id`: Gene ID of the corresponding gene
- `gene_name`: Gene name of the corresponding gene
- `gene_type`: Gene type of the corresponding gene

To get a guick summary of the eccDNA to see which genes are present, you can run the following command:

```bash
cut -f7 data/SRR413969.eccdna.annotated.tsv | sort | uniq -c | sort -nr
```

This provides counts on the number of eccDNA associated with each gene type.
Then you can extract each gene type and save it to a separate file as follows:

```bash
awk -F'\t' '$7 == "protein_coding" {print $6}' data/SRR413969.eccdna.annotated.tsv | sort | uniq > ecc_protein_genes.txt
```

This will create a file called `ecc_protein_genes.txt` with the list of protein-coding genes associated with eccDNA.
You can repeat this for any other gene type you are interested in.