#!/bin/bash
# get_ref_genome.sh
# Usage: bash get_ref_genome.sh

# Script to download and unzip the human genome reference (GRCh37)
# Author: Jason Hunter

# Download the full GRCh37 genome (GCF_000001405.13_GRCh37.p13)
echo "Downloading GRCh37 full genome reference (GCF_000001405.13)..."

curl -L -o data/human_genome/whole/GCF_000001405.13_GRCh37_genomic.fna.gz \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz"

# Check if download succeeded
if [ ! -f data/human_genome/whole/GCF_000001405.13_GRCh37_genomic.fna.gz ]; then
  echo "Error: Download failed."
  exit 1
fi


echo "Download Complete! Unzipping the genome file..."
# Unzip the downloaded file
gunzip data/human_genome/whole/GCF_000001405.13_GRCh37_genomic.fna.gz
# Check if unzip succeeded
if [ ! -f data/human_genome/whole/GCF_000001405.13_GRCh37_genomic.fna ]; then
  echo "Error: Unzipping failed."
  exit 1
fi
echo "Unzipping Complete! The genome file is ready for use."

