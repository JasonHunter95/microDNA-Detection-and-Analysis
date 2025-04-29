#!/bin/bash
# bwa_align_and_sort_by_chr.sh
# Usage: bash bwa_align_and_sort_by_chr.sh <chr_name> <sample_prefix>

# Script to align reads to a reference genome using BWA and sort the output BAM file
# Author: Jason Hunter

set -e

# Arguments
CHR=$1
SAMPLE=$2

# Paths
REF="data/human_genome/${CHR}/${CHR}.fna"
READ1="data/fastqfiles/${SAMPLE}_1.fastq"
READ2="data/fastqfiles/${SAMPLE}_2.fastq"
OUT_DIR="data/bams/${CHR}"
OUT_BAM="${OUT_DIR}/${SAMPLE}_${CHR}.sorted.bam"

# Check inputs
if [ -z "$CHR" ] || [ -z "$SAMPLE" ]; then
    echo "Usage: bash bwa_align_and_sort_by_chr.sh <chr_name> <sample_prefix>"
    echo "Example: bash bwa_align_and_sort_by_chr.sh chr1 SRR413984"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUT_DIR"

# Run BWA MEM and SAMTOOLS SORT
echo "Aligning ${SAMPLE} to ${CHR}..."
bwa mem -t 6 "$REF" "$READ1" "$READ2" | \
    samtools sort -@ 6 -m 2G -o "$OUT_BAM"

echo "Finished! Output BAM: $OUT_BAM"
