#!/bin/bash
# index_all_chr.sh
# Usage: bash index_all_chr.sh

# Script to unzip, index, and optionally rezip chromosome FASTA files
# Author: Jason Hunter

# Set the base directory
BASE_DIR="data/human_genome"

echo "Starting chromosome indexing pipeline..."
echo

# Step 1: Unzip all .fna.gz files if needed
echo "Unzipping all compressed FASTA files..."
find "$BASE_DIR" -name "*.fna.gz" -exec gunzip -k {} \;
echo "Unzipping complete."
echo

# Step 2: Index each .fna file
echo "Indexing each chromosome..."
for fasta in "$BASE_DIR"/chr*/chr*.fna; do
    if [ -f "$fasta" ]; then
        echo "Indexing $(basename "$fasta")..."
        bwa-mem2 index "$fasta"
    else
        echo "No .fna files found in $BASE_DIR/chr*/"
        exit 1
    fi
done
echo "Indexing complete."
echo

echo
echo "All done!"