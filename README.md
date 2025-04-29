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

Raw sequencing data can be fetched from the SRA database. The example dataset is from the NCBI SRA database, specifically the SRR413984 accession number.

```bash
fasterq-dump SRR413969 --split-files -O data/fastqfiles/
```

This will download the paired-end reads and save them in the `data/fastqfiles/` directory. The files will be named `SRR413984_1.fastq` and `SRR413984_2.fastq`.

The example dataset is a human sample, and the reference genome used in this example is the human genome (GRCh37). The reference genome file can be downloaded to the `data/human_genome/whole` directory as `GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna` via executing the following shell script:

Please note that this may take some time to download, as the reference genome is quite large.

```bash
bash src/utils/shell_scripts/get_ref_genome.sh
```

You can then use the following script to extract individual chromosomes from the genome downloaded using samtools faidx. This will create a folder for each chromosome in the genome, and subsequently saves an index file for each.
Once again, this may take a few minutes to run.

```bash
bash src/utils/shell_scripts/extract_chr.sh
```

## Alignment

The first step in the pipeline is to align the paired-end reads to the reference genome. This is done in two steps, first by indexing the reference genome and then aligning the reads to the indexed reference genome.

```bash
bash src/utils/shell_scripts/index_all_chr.sh
```

This will create an index file for each chromosome in the reference genome. The index files are used by BWA to align the reads to the reference genome.

It might also be useful to index the reference genome so that individual chromosomes can be accessed quickly.
A script which does this is also provided in the `src/utils/shell_scripts/` directory. It can be run as follows:

```bash
bash src/utils/shell_scripts/index_ref_genome.sh
```

This will create an index file for the reference genome, which is used by BWA to align the reads to the reference genome.
Alternatively, you can run the following command to index the reference genome using Samtools, which is a more common method:

```bash
samtools faidx data/human_genome/whole/GCF_000001405.13_GRCh37_genomic.fna
```

This can be done for each chromosome as well with the following command:

```bash
for i in {1..22} X Y M; do
data/human_genome/chr1
    samtools faidx data/human_genome/chr$i/chr$i.fna
done
```

Align the reads to the reference genome using BWA, and then sort the output BAM file using Samtools. The file is sent to the data/ folder. The `-t` option specifies the threadcount to use for BWA. I only used 4 threads for this because it was done locally on my mac. This can be increased if you are running on a more powerful machine.

```bash
bwa mem -t 6 data/human_genome/chr1/chr1.fna \
    data/fastqfiles/SRR413984_1.fastq data/fastqfiles/SRR413984_2.fastq | \
    samtools sort -@ 6 -m 2G -o data/bams/chr1/SRR413984_chr1.sorted.bam
```

If you want to run this for an individual chromosome, the following script can be used:

```bash
bash src/utils/shell_scripts/bwa_align_and_sort_by_chr.sh <chr_name> <sample_prefix>
```

Now, index the sorted BAM file (here just chr1) using Samtools:

```bash
samtools index data/bams/chr1/SRR413984_chr1.sorted.bam
```

Next, in order to prepare the BAM file for Circle-Map, we need to query sort the BAM file. This is done using the `samtools sort -n` command. The `-n` option sorts the reads by name, which is required for Circle-Map.

```bash
samtools sort -n -o data/bams/chr1/SRR413984_chr1.querysorted.bam data/bams/chr1/SRR413984_chr1.sorted.bam
```

Or for individual chromosomes, you can run the following script:

```bash
bash src/utils/shell_scripts/querysort_bam_by_chr.sh <chr_name> <sample_prefix>
```

## eccDNA Detection

The next step is to detect eccDNA using the Circle-Map tool. The Circle-Map pipeline consists of two main steps: ReadExtractor and Realign.
The ReadExtractor step extracts the reads that are most statistically likely to be circular, and the Realign step re-aligns the extracted reads to the reference genome.

```bash
Circle-Map ReadExtractor -i data/bams/chr1/SRR413984_chr1.querysorted.bam -o data/bams/chr1/SRR413984_chr1.candidates.bam
```

Then sort and index the candidates.bam file using Samtools:

```bash
samtools sort -o data/bams/chr1/SRR413984_chr1.candidates.sorted.bam data/bams/chr1/SRR413984_chr1.candidates.bam
samtools index data/bams/chr1/SRR413984_chr1.candidates.sorted.bam
```

Heres a script to do this for individual chromosomes:

```bash
bash src/utils/shell_scripts/sort_and_index_candidates_by_chr.sh chr< replace_with_your_chr_number > < sample_id >
```

Now we can run the Realign step. This step takes the queryname sorted BAM file and the original sorted BAM file as input, and outputs a BED file with the eccDNA coordinates.

```bash
Circle-Map Realign \
  -i data/bams/chr1/SRR413984_chr1.candidates.sorted.bam \
  -qbam data/bams/chr1/SRR413984_chr1.querysorted.bam \
  -sbam data/bams/chr1/SRR413984_chr1.sorted.bam \
  -fasta data/human_genome/chr1/chr1.fna \
  -o data/beds/chr1/SRR413984_chr1.eccdna.bed \
  --split 2 \
  --threads 8
```

## Post-processing

Run the following script to clean the output files:

```bash
python3 src/utils/clean_bed.py \
  -i data/beds/chr1/SRR413984_chr1.eccdna.bed \
  -o data/beds/chr1/SRR413984_chr1.eccdna.cleaned.bed
```

## Annotation

  1. Annotate the eccDNA using the Gencode annotation file.

      ```bash
        curl -L -o data/gtfs/human_genome/gencode.v19.annotation.gtf.gz \
           https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
        gunzip data/gtfs/human_genome/gencode.v19.annotation.gtf.gz
      ```

  2. Convert the GTF file to BED format using the `gtf2bed` command. This creates a BED file with the Gencode annotation.

      ```bash
        gtf2bed < data/gtfs/human_genome/gencode.v19.annotation.gtf > data/beds/whole_genomes/human_genome/gencode.v19.annotation.bed
      ```

      This will create a BED file with the Gencode annotation in the `data/beds/whole_genomes/human_genome/` directory.
      To extract the chromosome 1 annotation, you can use the following command:

      ```bash
          awk -F'\t' '$1 == "chr1"' data/beds/whole_genomes/human_genome/gencode.v19.annotation.bed > data/beds/chr1/gencode_chr1.v19.annotation.bed
      ```

      This will create a BED file with the Gencode annotation for chromosome 1 in the `data/beds/chr1/` directory.
      Then you can run the following command to filter only "gene" features from the Gencode annotation:

      ```bash
          awk '$8 == "gene"' data/beds/chr1/gencode_chr1.v19.annotation.bed > data/beds/chr1/gencode_chr1_only_genes.bed
      ```

      We also need to sort our data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.bed to give closest the proper inputs:

      ```bash
        bedtools sort -i data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.bed > data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.sorted.bed
      ```

  3. Convert the eccDNA BED file to UCSC format. This is necessary for the `bedtools intersect` command to work. The Gencode annotation uses `chr1`, `chr2`, etc. as chromosome names, while the CircleMap output uses `NC_000001.10`, `NC_000002.11`, etc. as chromosome names.

      ```bash
        sed 's/^NC_000001.10/chr1/' data/beds/chr1/SRR413984_chr1.eccdna.cleaned.bed > data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.bed
      ```

  4. Intersect the eccDNA BED file with the Gencode annotation BED file using `bedtools intersect`. This will create a new BED file with the eccDNA coordinates and their corresponding gene annotations.

      ```bash
        bedtools closest \
          -a data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.sorted.bed \
          -b data/beds/chr1/gencode_chr1_only_genes.bed \       
          -d > data/beds/chr1/SRR413984_chr1.eccdna.closest_genes.bed
      ```

  5. Generate a clean output of the eccDNA with gene annotations using the `annotate_eccdna_genes.py` script. This will create a new TSV file with the eccDNA coordinates and their corresponding gene annotations.

      ```bash
        python3 src/utils/annotate_eccdna_closest.py
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
cut -f7 data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv | sort | uniq -c | sort -nr
```

This provides counts on the number of eccDNA associated with each gene type.
Then you can extract each gene type and save it to a separate file as follows:

```bash
awk '$7 == "protein_coding"' data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv > data/results/chr1/SRR413984_chr1_protein_coding_genes.tsv
awk '$7 == "lncRNA"' data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv > data/results/chr1/SRR413984_chr1_lncRNA_genes.tsv
awk '$7 == "miRNA"' data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv > data/results/chr1/SRR413984_chr1_miRNA_genes.tsv
awk '$7 == "snRNA"' data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv > data/results/chr1/SRR413984_chr1_snRNA_genes.tsv
awk '$7 == "snoRNA"' data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv > data/results/chr1/SRR413984_chr1_snoRNA_genes.tsv
awk '$7 == "rRNA"' data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv > data/results/chr1/SRR413984_chr1_rRNA_genes.tsv
awk '$7 == "pseudogene"' data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv > data/results/chr1/SRR413984_chr1_pseudogene_genes.tsv
awk '$7 == "miRNA"' data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv > data/results/chr1/SRR413984_chr1_miRNA_genes.tsv
awk '$7 == "tRNA"' data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv > data/results/chr1/SRR413984_chr1_tRNA_genes.tsv
```
