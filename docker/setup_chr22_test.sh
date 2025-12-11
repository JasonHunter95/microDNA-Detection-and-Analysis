#!/bin/bash
# setup_chr22_test.sh - Set up chromosome 22 subset for low-memory local testing of the CircleMap pipeline
#
# For systems with <32GB RAM that cannot build the full genome BWA-MEM2 index.
# Uses chr22 (smallest autosome) as a representative subset for end-to-end testing.
#
# Prerequisites:
#   - Full reference genome already downloaded (run setup_data.sh --skip-index first)
#   - Annotations already prepared
#   - Docker image built
#
# Usage:
#   ./docker/setup_chr22_test.sh
#
# This script will:
#   1. Extract chr22 from the full reference genome
#   2. Create a BWA-MEM2 index for chr22 (fits in 16GB RAM)
#   3. Create chr22-specific annotation subset

set -e

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

IMAGE_NAME="microdna-circlemap:latest"
FULL_REFERENCE="data/inputs/references/genome/GCF_000001405.13_GRCh37_genomic.fna"
CHR22_REFERENCE="data/inputs/references/genome/chr22.fna"
ANNOTATIONS_DIR="data/inputs/references/annotations"
FULL_GENE_BED="${ANNOTATIONS_DIR}/gencode.v19.annotation.genes.sorted.bed"
CHR22_GENE_BED="${ANNOTATIONS_DIR}/gencode.v19.chr22.genes.sorted.bed"

# NCBI accession for chr22 in GRCh37
CHR22_ACCESSION="NC_000022.10"

echo "========================================"
echo "  Chr22 Subset Setup (Low-Memory Mode)"
echo "========================================"
echo ""

# Validate prerequisites
if [ ! -f "$FULL_REFERENCE" ]; then
    echo "ERROR: Full reference genome not found."
    echo "       Run './docker/setup_data.sh --skip-index' first to download it."
    exit 1
fi

if [ ! -f "$FULL_GENE_BED" ]; then
    echo "ERROR: Annotations not found."
    echo "       Run './docker/setup_data.sh --skip-reference --skip-index' first."
    exit 1
fi

# Build Docker image if needed
if ! docker image inspect "$IMAGE_NAME" &> /dev/null; then
    echo "Building Docker image..."
    docker build -t "$IMAGE_NAME" -f "$PROJECT_DIR/docker/Dockerfile.circlemap" "$PROJECT_DIR/docker/"
fi

# Helper to run commands in Docker
run_docker() {
    docker run --rm -v "$PROJECT_DIR:/workspace" -w /workspace "$IMAGE_NAME" "$@"
}

# Step 1: Extract chr22 from reference
if [ -f "$CHR22_REFERENCE" ]; then
    echo "[1/3] Chr22 reference already exists, skipping extraction..."
else
    echo "[1/3] Extracting chr22 from full reference..."
    
    # Use samtools faidx to extract chr22
    # First ensure faidx exists
    if [ ! -f "${FULL_REFERENCE}.fai" ]; then
        echo "      Creating FASTA index..."
        run_docker samtools faidx "$FULL_REFERENCE"
    fi
    
    # Extract chr22 using the NCBI accession
    run_docker samtools faidx "$FULL_REFERENCE" "$CHR22_ACCESSION" > "$CHR22_REFERENCE"
    
    # Rename header to simpler "chr22" for compatibility
    sed -i.bak "s/>${CHR22_ACCESSION}.*/>chr22/" "$CHR22_REFERENCE" && rm -f "${CHR22_REFERENCE}.bak"
    
    echo "      ✓ Chr22 extracted ($(du -h "$CHR22_REFERENCE" | cut -f1))"
fi

# Step 2: Create BWA-MEM2 index for chr22
if [ -f "${CHR22_REFERENCE}.bwt.2bit.64" ]; then
    echo "[2/3] Chr22 BWA-MEM2 index already exists, skipping..."
else
    echo "[2/3] Creating BWA-MEM2 index for chr22 (this takes ~5-10 minutes)..."
    run_docker bwa-mem2 index "$CHR22_REFERENCE"
    run_docker samtools faidx "$CHR22_REFERENCE"
    echo "      ✓ BWA-MEM2 index created"
fi

# Step 3: Create chr22-specific annotations
if [ -f "$CHR22_GENE_BED" ]; then
    echo "[3/3] Chr22 annotations already exist, skipping..."
else
    echo "[3/3] Extracting chr22-specific annotations..."
    
    # Filter annotations for chr22 (both "chr22" and NCBI accession formats)
    grep -E "^(chr22|${CHR22_ACCESSION})\s" "$FULL_GENE_BED" > "$CHR22_GENE_BED" || true
    
    # If using GENCODE which uses chr22 naming, also normalize
    if [ ! -s "$CHR22_GENE_BED" ]; then
        echo "      Warning: No chr22 annotations found. The annotation file may use different naming."
        echo "      Attempting alternate extraction..."
        awk '$1 == "chr22" || $1 == "'$CHR22_ACCESSION'"' "$FULL_GENE_BED" > "$CHR22_GENE_BED" || true
    fi
    
    COUNT=$(wc -l < "$CHR22_GENE_BED" | tr -d ' ')
    echo "      ✓ Extracted $COUNT gene annotations for chr22"
fi

echo ""
echo "========================================"
echo "  ✅ Chr22 subset setup complete!"
echo "========================================"
echo ""
echo "Files created:"
echo "  - $CHR22_REFERENCE"
echo "  - ${CHR22_REFERENCE}.bwt.2bit.64 (BWA-MEM2 index)"
echo "  - $CHR22_GENE_BED"
echo ""
echo "Next steps:"
echo "  1. Align reads to chr22:"
echo "     ./docker/run_alignment.sh SRR413984 4 chr22"
echo ""
echo "  2. Run full pipeline:"
echo "     ./docker/run_full_pipeline.sh \\"
echo "         data/intermediate/SRR413984.chr22.sorted.bam \\"
echo "         $CHR22_REFERENCE \\"
echo "         $CHR22_GENE_BED \\"
echo "         data/outputs/SRR413984_chr22"
