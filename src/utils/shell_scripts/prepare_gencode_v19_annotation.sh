#!/bin/bash
# Filename: prepare_gencode_v19_annotation.sh
# Purpose: Download and prepare Gencode v19 annotation for eccDNA pipeline

set -e # Exit immediately if a command exits with a non-zero status

# --- Directories ---
ANNOTATIONS_DIR="data/inputs/references/annotations"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
GTF_GZ_FILE="${ANNOTATIONS_DIR}/gencode.v19.annotation.gtf.gz"
GTF_FILE="${ANNOTATIONS_DIR}/gencode.v19.annotation.gtf"
FULL_BED_FILE="${ANNOTATIONS_DIR}/gencode.v19.annotation.bed"
GENE_BED_FILE="${ANNOTATIONS_DIR}/gencode.v19.annotation.genes.bed"
# This is the final output file needed by the main pipeline
SORTED_GENE_BED_FILE="${ANNOTATIONS_DIR}/gencode.v19.annotation.genes.sorted.bed"

echo "Creating directories..."
mkdir -p "${ANNOTATIONS_DIR}"

# --- Check if final file already exists ---
if [ -f "${SORTED_GENE_BED_FILE}" ]; then
    echo "✅ Final sorted gene BED file already exists: ${SORTED_GENE_BED_FILE}"
    echo "Skipping preparation."
    exit 0
fi

echo "--- Starting Annotation Preparation ---"

# 1. Download (only if GTF doesn't exist)
if [ ! -f "${GTF_FILE}" ] && [ ! -f "${GTF_GZ_FILE}" ]; then
    echo "Downloading GTF..."
    curl -L -o "${GTF_GZ_FILE}" "${GTF_URL}"
else
    echo "GTF file or gzipped GTF found, skipping download."
fi

# 2. Unzip (only if unzipped version doesn't exist)
if [ ! -f "${GTF_FILE}" ] && [ -f "${GTF_GZ_FILE}" ]; then
    echo "Unzipping GTF..."
    gunzip -f "${GTF_GZ_FILE}" # Use -f to overwrite just in case
else
    echo "Unzipped GTF found or GZ file missing, skipping unzip."
fi

# Ensure GTF file exists before proceeding
if [ ! -f "${GTF_FILE}" ]; then
    echo "Error: GTF file not found at ${GTF_FILE} after download/unzip steps." >&2
    exit 1
fi

# 3. Convert GTF to BED
echo "Converting GTF to BED..."
gtf2bed < "${GTF_FILE}" > "${FULL_BED_FILE}"

# 4. Filter for Genes (Check column $8 for your gtf2bed output)
echo "Filtering for gene features..."
awk -F'\t' '$8 == "gene"' "${FULL_BED_FILE}" > "${GENE_BED_FILE}"

# 5. Sort Gene BED
echo "Sorting gene BED file..."
bedtools sort -i "${GENE_BED_FILE}" > "${SORTED_GENE_BED_FILE}"

echo "--- Annotation preparation complete ---"
echo "✅ Final file: ${SORTED_GENE_BED_FILE}"

# Optional: Clean up intermediate files
# echo "Cleaning up intermediate files..."
# rm "${FULL_BED_FILE}" "${GENE_BED_FILE}"

exit 0