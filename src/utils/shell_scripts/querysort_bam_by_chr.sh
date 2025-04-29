#!/bin/bash
# querysort_bam_by_chr.sh
# Usage: bash querysort_bam_by_chr.sh <chr_name> <sample_prefix>

# Script to queryname sort BAM files by chromosome
# Author: Jason Hunter

set -e

# Arguments
CHR=$1
SAMPLE=$2

# Paths
INPUT_BAM="data/bams/${CHR}/${SAMPLE}_${CHR}.sorted.bam"
OUTPUT_BAM="data/bams/${CHR}/${SAMPLE}_${CHR}.querysorted.bam"

# Check inputs
if [ -z "$CHR" ] || [ -z "$SAMPLE" ]; then
    echo "Usage: bash querysort_bam_by_chr.sh <chr_name> <sample_prefix>"
    echo "Example: bash querysort_bam_by_chr.sh chr1 SRR413984"
    exit 1
fi

# Check if input BAM exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file not found: $INPUT_BAM"
    exit 1
fi

# Queryname sort
echo "Sorting $INPUT_BAM by query name..."
samtools sort -n -@ 6 -m 2G -o "$OUTPUT_BAM" "$INPUT_BAM"

echo "Finished! Output BAM: $OUTPUT_BAM"
