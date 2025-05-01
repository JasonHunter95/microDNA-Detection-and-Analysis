# microDNA-Detection-and-Analysis for Windows 11

This guide provides step-by-step instructions for setting up and running the microDNA-Detection-and-Analysis pipeline on Windows 11.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Installation](#installation)
3. [Data Acquisition](#data-acquisition)
4. [Alignment and Indexing](#alignment-and-indexing)
5. [eccDNA Detection](#eccdna-detection)
6. [Post-processing](#post-processing)
7. [Annotation](#annotation)
8. [Results](#results)
9. [Troubleshooting](#troubleshooting)

## Prerequisites

For Windows 11, we'll need to set up the following:

- **Windows Subsystem for Linux (WSL2)** - Provides a Linux environment on Windows
- **Ubuntu** - Our Linux distribution of choice for WSL
- **Conda** (Miniconda) - Package manager for Python environments
- All other bioinformatics tools will be installed through conda

## Installation

### 1. Install WSL2 and Ubuntu

Open PowerShell as Administrator and run:

```powershell
wsl --install
```

This command installs WSL2 with Ubuntu by default. After installation completes, restart your computer.

After restart, Ubuntu will open automatically and ask you to create a username and password. Follow the prompts to complete setup.

### 2. Set up conda in WSL Ubuntu

Launch Ubuntu from the Start menu and run the following commands:

```bash
# Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh

# Make the installer executable
chmod +x miniconda.sh

# Run the installer
./miniconda.sh

# Follow the prompts to install Miniconda
# Say yes when asked to initialize Miniconda

# Reload your bash profile
source ~/.bashrc

# Verify conda installation
conda --version
```

### 3. Clone the Repository

In your Ubuntu WSL terminal:

```bash
# Navigate to your Windows Desktop (accessible through /mnt/c/Users/...)
cd /mnt/c/Users/Jason\ Hunter/Desktop

# Clone the repository
git clone https://github.com/JasonHunter95/microDNA-Detection-and-Analysis.git
cd microDNA-Detection-and-Analysis
```

### 4. Create and Activate the conda Environment

```bash
# Create environment from the provided yml file
conda env create -f environment.yml

# Activate the environment
conda activate circlemap-env
```

If the environment.yml file doesn't work properly on Windows/WSL, create the environment manually:

```bash
# Create a new conda environment for Circle-Map
conda create -n circlemap-env python=3.6
conda activate circlemap-env

# Install required packages
conda install -c bioconda circle-map
conda install -c bioconda bwa-mem2
conda install -c bioconda samtools
conda install -c bioconda bedtools
conda install -c bioconda sra-tools
conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed

# Install additional Python packages if needed
pip install pandas numpy
```

### 5. Create Directory Structure

```bash
mkdir -p data/fastqfiles
mkdir -p data/human_genome/whole
mkdir -p data/human_genome/chr1
mkdir -p data/bams/chr1
mkdir -p data/beds/chr1
mkdir -p data/gtfs/human_genome
mkdir -p data/beds/whole_genomes/human_genome
mkdir -p data/results/chr1
mkdir -p src/utils
```

## Data Acquisition

### 1. Download Reference Data

```bash
# Download human reference genome (chromosome 1)
cd /mnt/c/Users/Jason\ Hunter/microDNA-Detection-and-Analysis
mkdir -p data/human_genome/chr1

# Download chromosome 1 of the human genome
wget -O data/human_genome/chr1/chr1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr1.fna.gz
gunzip data/human_genome/chr1/chr1.fna.gz

# Download the Gencode annotation file
mkdir -p data/gtfs/human_genome
wget -O data/gtfs/human_genome/gencode.v19.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip data/gtfs/human_genome/gencode.v19.annotation.gtf.gz
```

### 2. Download Example Data

```bash
# Download example data (this will take some time)
cd /mnt/c/Users/Jason\ Hunter/microDNA-Detection-and-Analysis
mkdir -p data/fastqfiles
fasterq-dump SRR413984 --split-files -O data/fastqfiles/
```

### 3. Prepare Annotation Files

```bash
# Convert GTF to BED format
mkdir -p data/beds/whole_genomes/human_genome

# Extract only chromosome 1 from the GTF file
grep -P "^chr1\t" data/gtfs/human_genome/gencode.v19.annotation.gtf > data/gtfs/human_genome/gencode.v19.chr1.gtf

# Convert chromosome 1 GTF to BED format
mkdir -p data/beds/chr1

# Install bedops
conda install -c bioconda bedops
convert2bed -i gtf < data/gtfs/human_genome/gencode.v19.chr1.gtf > data/beds/chr1/gencode_chr1_only_genes.bed

# Sort the BED file
bedtools sort -i data/beds/chr1/gencode_chr1_only_genes.bed > data/beds/chr1/gencode_chr1_only_genes.sorted.bed
```

## Alignment and Indexing

```bash
# Index the chromosome 1 reference
bwa-mem2 index data/human_genome/chr1/chr1.fna

# Step 1: Align reads and save to SAM file
bwa-mem2 mem -t 12 data/human_genome/chr1/chr1.fna \
    data/fastqfiles/SRR413984_1.fastq data/fastqfiles/SRR413984_2.fastq \
    > data/bams/chr1/SRR413984_chr1.sam

# Step 2: Convert SAM to sorted BAM
samtools sort -@ 12 -m 1G \
    -o data/bams/chr1/SRR413984_chr1.sorted.bam \
    data/bams/chr1/SRR413984_chr1.sam

# Optional: Remove the intermediate SAM file to save space
rm data/bams/chr1/SRR413984_chr1.sam

# Index the sorted BAM file
samtools index data/bams/chr1/SRR413984_chr1.sorted.bam

# Query-sort the BAM file for Circle-Map
samtools sort -n -o data/bams/chr1/SRR413984_chr1.querysorted.bam data/bams/chr1/SRR413984_chr1.sorted.bam

# Try to index for Circle-Map but I think it will fail (is there a way to do this?)
# Circle-Map seems to expect a query-sorted BAM file but you can't index one? LIKE WTF
samtools index data/bams/chr1/SRR413984_chr1.querysorted.bam
```

## eccDNA Detection

```bash
# Extract candidate circular reads
Circle-Map ReadExtractor -i data/bams/chr1/SRR413984_chr1.querysorted.bam -o data/bams/chr1/SRR413984_chr1.candidates.bam

# Sort and index candidate reads
samtools sort -o data/bams/chr1/SRR413984_chr1.candidates.sorted.bam data/bams/chr1/SRR413984_chr1.candidates.bam
samtools index data/bams/chr1/SRR413984_chr1.candidates.sorted.bam

# Run Circle-Map Realign to detect eccDNA
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

To clean up the eccDNA output, create a Python script for cleaning:

```bash
# Create a clean_bed.py script
cat > src/utils/clean_bed.py << 'EOL'
#!/usr/bin/env python3

import argparse
import pandas as pd

def clean_bed(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'circ_reads', 'disc_reads', 'coverage'])
    
    # Filter by size (100-400 bp)
    df['size'] = df['end'] - df['start']
    df_filtered = df[(df['size'] >= 100) & (df['size'] <= 400)]
    
    # Save to file
    df_filtered.to_csv(output_file, sep='\t', index=False, header=False)
    print(f"Filtered {len(df) - len(df_filtered)} entries. Remaining: {len(df_filtered)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Clean and filter eccDNA BED file')
    parser.add_argument('-i', '--input', required=True, help='Input BED file')
    parser.add_argument('-o', '--output', required=True, help='Output BED file')
    args = parser.parse_args()
    
    clean_bed(args.input, args.output)
EOL

# Make the script executable
chmod +x src/utils/clean_bed.py

# Run the cleaning script
python3 src/utils/clean_bed.py \
  -i data/beds/chr1/SRR413984_chr1.eccdna.bed \
  -o data/beds/chr1/SRR413984_chr1.eccdna.cleaned.bed
```

## Annotation

```bash
# Convert the eccDNA BED file to UCSC format
sed 's/^chr1/chr1/' data/beds/chr1/SRR413984_chr1.eccdna.cleaned.bed > data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.bed

# Sort the eccDNA BED file
bedtools sort -i data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.bed > data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.sorted.bed

# Find closest genes to eccDNA
bedtools closest \
  -a data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.sorted.bed \
  -b data/beds/chr1/gencode_chr1_only_genes.sorted.bed \
  -d > data/beds/chr1/SRR413984_chr1.eccdna.closest_genes.bed
```

Create an annotation script:

```bash
# Create a script for annotating eccDNA
cat > src/utils/annotate_eccdna_closest.py << 'EOL'
#!/usr/bin/env python3

import pandas as pd

# Load eccDNA data
eccdna_file = "data/beds/chr1/SRR413984_chr1.eccdna.closest_genes.bed"
eccdna_df = pd.read_csv(eccdna_file, sep='\t', header=None,
                        names=['ecc_chrom', 'ecc_start', 'ecc_end', 'ecc_name', 'ecc_score', 'ecc_strand',
                               'circ_reads', 'disc_reads', 'coverage', 
                               'gene_chrom', 'gene_start', 'gene_end', 'gene_name', 'gene_score', 'gene_strand',
                               'distance'])

# Clean and prepare output
result_df = pd.DataFrame({
    'ecc_id': eccdna_df['ecc_name'],
    'chrom': eccdna_df['ecc_chrom'],
    'start': eccdna_df['ecc_start'],
    'end': eccdna_df['ecc_end'],
    'size': eccdna_df['ecc_end'] - eccdna_df['ecc_start'],
    'gene_name': eccdna_df['gene_name'],
    'distance': eccdna_df['distance']
})

# Save to file
output_file = "data/results/chr1/SRR413984_chr1_annotated.eccdna.tsv"
result_df.to_csv(output_file, sep='\t', index=False)
print(f"Saved annotated eccDNA to {output_file}")
EOL

# Make the script executable
chmod +x src/utils/annotate_eccdna_closest.py

# Run the annotation script
python3 src/utils/annotate_eccdna_closest.py
```

## Results

Generate summary statistics:

```bash
# Count eccDNA by gene
cut -f6 data/results/chr1/SRR413984_chr1_annotated.eccdna.tsv | sort | uniq -c | sort -nr > data/results/chr1/SRR413984_chr1_gene_counts.txt

# Generate size distribution stats
python -c "
import pandas as pd
df = pd.read_csv('data/results/chr1/SRR413984_chr1_annotated.eccdna.tsv', sep='\t')
print('Size statistics:')
print(df['size'].describe())
"
```

## Troubleshooting

### Common Issues and Solutions

1. **WSL File Permission Issues**
   - When accessing files in the Windows file system from WSL, you might encounter permission issues
   - Solution: Work directly within the WSL file system or use `chmod` to set appropriate permissions

2. **Conda Environment Activation**
   - If you get `CommandNotFoundError: Your shell has not been properly configured to use 'conda activate'`
   - Solution: Run `conda init bash` and restart your terminal

3. **Running Out of Memory**
   - For large genomes, reduce the number of threads or increase WSL memory allocation
   - Edit `.wslconfig` file in your Windows user directory to increase memory: `memory=8GB`

4. **Improving WSL Performance**
   - Move your working directory inside the WSL filesystem for better performance
   - Example: `mkdir ~/microDNA-project && cp -r /mnt/c/Users/Jason\ Hunter/Desktop/microDNA-Detection-and-Analysis/* ~/microDNA-project/`

5. **Installing BEDOPS Tools**
   - If `gtftogenepred` or other tools aren't available: `conda install -c bioconda ucsc-gtftogenepred ucsc-genepredtobed bedops`

For additional help, consult the original repository documentation or create an issue on GitHub.