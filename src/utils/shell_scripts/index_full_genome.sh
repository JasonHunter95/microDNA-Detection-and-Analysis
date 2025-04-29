#!/bin/bash
# index_full_genome.sh
# Usage: bash index_full_genome.sh

# Script to index the entire genome FASTA file
# Author: Jason Hunter

# Set the path to the full genome FASTA
GENOME_FASTA="data/human_genome/whole/GCF_000001405.13_GRCh37_genomic.fna"

echo "Starting full genome indexing..."
echo

# Step 1: Check if the genome FASTA exists
if [ ! -f "$GENOME_FASTA" ]; then
    echo "Error: Genome FASTA file not found at $GENOME_FASTA"
    exit 1
fi

# Step 2: Index the genome
echo "Indexing $(basename "$GENOME_FASTA")..."
bwa index "$GENOME_FASTA"

echo
echo "Genome indexing complete!"
