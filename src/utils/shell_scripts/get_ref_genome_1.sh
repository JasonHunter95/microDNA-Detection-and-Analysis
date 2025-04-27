#!/bin/bash

# Download only chromosome 1 (GRCh37, NC_000001.10)
echo "Downloading GRCh37 chr1 (NC_000001.10) reference..."

curl -L -o data/GRCh37_chr1.fna.gz \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.fna.gz"

# check if download succeeded
if [ ! -f data/GRCh37_chr1.fna.gz ]; then
  echo "Error: Download failed."
  exit 1
fi

# extract it
echo "Extracting reference..."
gunzip -c data/GRCh37_chr1.fna.gz > data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna

echo "All done!"