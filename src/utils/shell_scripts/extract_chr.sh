#!/bin/bash

# Usage: bash src/utils/shell_scripts/extract_chr.sh

INPUT_FASTA="data/human_genome/GCF_000001405.13_GRCh37_genomic.fna"
BASE_OUTPUT_DIR="data/human_genome"

# Index if missing
if [ ! -f "${INPUT_FASTA}.fai" ]; then
  echo "Indexing FASTA file..."
  samtools faidx "$INPUT_FASTA"
fi

# Extract list of primary chromosomes dynamically
grep "^>" "$INPUT_FASTA" | while read line; do
  accession=$(echo "$line" | cut -d' ' -f1 | sed 's/>//')
  desc=$(echo "$line" | cut -d' ' -f2-)

  # ⬇️ Check accession prefix: only NC_... or NC_012920.1
  if [[ ! "$accession" =~ ^NC_ ]]; then
    echo "Skipping non-primary scaffold: $accession"
    continue
  fi

  # Detect if it's a primary chromosome or mitochondrion
  if [[ "$desc" =~ chromosome\ ([0-9]+) ]]; then
    chr="chr${BASH_REMATCH[1]}"
  elif [[ "$desc" =~ chromosome\ (X) ]]; then
    chr="chrX"
  elif [[ "$desc" =~ chromosome\ (Y) ]]; then
    chr="chrY"
  elif [[ "$desc" =~ mitochondrion ]]; then
    chr="chrM"
  else
    echo "Skipping unrelated sequence: $accession"
    continue
  fi

  echo "Extracting $chr from $accession..."

  chr_dir="${BASE_OUTPUT_DIR}/${chr}"
  chr_file="${chr_dir}/${chr}.fna.gz"

  mkdir -p "$chr_dir"

  if samtools faidx "$INPUT_FASTA" "$accession" | gzip > "$chr_file"; then
    echo "Saved $chr to $chr_file (compressed)"
  else
    echo "Warning: Failed to extract $accession."
  fi
done

echo "All primary chromosomes extracted and compressed!"
