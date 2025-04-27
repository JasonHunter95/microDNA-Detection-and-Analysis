#!/bin/bash

# download the reference genome (only chromosome 1 from GRCh37)
# this can be modified to download other chromosomes or assemblies in the future
echo "Downloading GRCh37 chr1 reference..."
wget -O data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37.p13/GCF_000001405.13_GRCh37.p13_genomic.fna.gz"

# extract it if it's compressed, optionally
echo "Extracting reference..."
gunzip -c data/GCF_000001405.13_GRCh37.p13_genomic.fna.gz > data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna

echo "Done!"