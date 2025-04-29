#!/bin/bash
# sort_and_index_candidates_by_chr.sh
# Usage: bash sort_and_index_candidates_by_chr.sh <chr_name> <sample_prefix>

# Script to sort and index candidate BAM files by chromosome
# Author: Jason Hunter

set -e

# Arguments
CHR=$1
SAMPLE=$2

# Paths
INPUT_BAM="data/bams/${CHR}/${SAMPLE}_${CHR}.candidates.bam"
SORTED_BAM="data/bams/${CHR}/${SAMPLE}_${CHR}.candidates.sorted.bam"

# Check inputs
if [ -z "$CHR" ] || [ -z "$SAMPLE" ]; then
    echo "Usage: bash sort_index_candidates_by_chr.sh <chr_name> <sample_prefix>"
    echo "Example: bash sort_index_candidates_by_chr.sh chr1 SRR413984"
    exit 1
fi

# Check if input BAM exists
if [ ! -f "$INPUT_BAM" ]; then
    echo "Error: Input BAM file not found: $INPUT_BAM"
    exit 1
fi

# Sort the candidates BAM
echo "Sorting $INPUT_BAM..."
samtools sort -@ 6 -m 2G -o "$SORTED_BAM" "$INPUT_BAM"

# Index the sorted BAM
echo "Indexing $SORTED_BAM..."
samtools index "$SORTED_BAM"

echo "Finished! Output BAM: $SORTED_BAM"
echo "Finished! Output BAI: ${SORTED_BAM}.bai"
