#!/bin/bash
# post_pipeline_cleanup.sh
# Usage: bash post_pipeline_cleanup.sh
# Safely removes intermediate files to reduce local storage use
# Author: Jason Hunter

set -e

echo "Starting post-processing cleanup..."

# 1. Delete very large intermediate BAMs
echo "Deleting intermediate BAM files..."
rm -v data/bams/chr1/SRR413984_chr1.querysorted.bam || true
rm -v data/bams/chr1/SRR413984_chr1.sorted.bam || true
rm -v data/bams/chr1/SRR413984_chr1.candidates.bam || true
rm -v data/bams/chr1/SRR413984_chr1.candidates.sorted.bam || true

# 2. Remove huge Circle-Map annotation BED
echo "Deleting large Circle-Map with_genes BED..."
rm -v data/beds/chr1/SRR413984_chr1.eccdna.with_genes.bed || true

# 3. Optionally compress FASTQ files (comment out if you want to keep raw)
echo "Compressing FASTQ files..."
gzip -v data/fastqfiles/*.fastq || true

# 4. Optionally remove full reference genome if working per-chromosome only
echo "Deleting full reference genome (keeping per-chromosome files)..."
rm -v data/human_genome/whole/GCF_000001405.13_GRCh37_genomic.fna* || true

# 5. Compress Gencode GTF file
echo "Compressing Gencode GTF annotation..."
gzip -v data/gtfs/gencode.v19.annotation.gtf || true

echo "Cleanup complete! Disk space should be greatly reduced."
