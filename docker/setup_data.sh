#!/bin/bash
# setup_data.sh - Download all required data for the microDNA pipeline
#
# All steps run in Docker for maximum reproducibility.
#
# Usage:
#   ./docker/setup_data.sh [options]
#
# Options:
#   --sample ID      Download FASTQ from SRA (e.g., SRR413984)
#   --skip-reference Skip reference genome download
#   --skip-annotation Skip annotation preparation
#   --skip-index     Skip BWA indexing
#
# Example:
#   ./docker/setup_data.sh --sample SRR413984

set -e

# Parse arguments
SAMPLE_ID=""
SKIP_REFERENCE=false
SKIP_ANNOTATION=false
SKIP_INDEX=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample) SAMPLE_ID="$2"; shift 2 ;;
        --skip-reference) SKIP_REFERENCE=true; shift ;;
        --skip-annotation) SKIP_ANNOTATION=true; shift ;;
        --skip-index) SKIP_INDEX=true; shift ;;
        -h|--help)
            echo "Usage: $0 [--sample ID] [--skip-reference] [--skip-annotation] [--skip-index]"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

IMAGE_NAME="microdna-circlemap:latest"
REFERENCE="data/inputs/references/genome/GCF_000001405.13_GRCh37_genomic.fna"
ANNOTATIONS_DIR="data/inputs/references/annotations"
SORTED_GENE_BED="${ANNOTATIONS_DIR}/gencode.v19.annotation.genes.sorted.bed"

echo "========================================"
echo "  microDNA Data Setup (Docker)"
echo "========================================"
echo ""

# Build Docker image if needed (uses existing Circle-Map image with all tools)
if ! docker image inspect "$IMAGE_NAME" &> /dev/null; then
    echo "Building Docker image..."
    docker build -t "$IMAGE_NAME" -f "$PROJECT_DIR/docker/Dockerfile.circlemap" "$PROJECT_DIR/docker/"
fi

# Helper to run commands in Docker
run_docker() {
    docker run --rm -v "$PROJECT_DIR:/workspace" -w /workspace "$IMAGE_NAME" "$@"
}

# Step 1: Reference genome
if [ "$SKIP_REFERENCE" = false ]; then
    if [ -f "$REFERENCE" ]; then
        echo "[1/4] Reference genome already exists, skipping..."
    else
        echo "[1/4] Downloading GRCh37 reference genome (~3GB)..."
        mkdir -p data/inputs/references/genome
        
        run_docker bash -c "
            curl -L -o ${REFERENCE}.gz \
                'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz' && \
            gunzip ${REFERENCE}.gz
        "
        echo "    ✓ Reference genome downloaded"
    fi
else
    echo "[1/4] Skipping reference genome (--skip-reference)"
fi

# Step 2: Annotations
if [ "$SKIP_ANNOTATION" = false ]; then
    if [ -f "$SORTED_GENE_BED" ]; then
        echo "[2/4] Annotations already exist, skipping..."
    else
        echo "[2/4] Preparing GENCODE v19 annotations..."
        mkdir -p "$ANNOTATIONS_DIR"
        
        GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
        GTF_FILE="${ANNOTATIONS_DIR}/gencode.v19.annotation.gtf"
        
        run_docker bash -c "
            curl -L -o ${GTF_FILE}.gz '$GTF_URL' && \
            gunzip -f ${GTF_FILE}.gz && \
            gtf2bed < '$GTF_FILE' > '${ANNOTATIONS_DIR}/gencode.v19.annotation.bed' && \
            awk -F'\t' '\$8 == \"gene\"' '${ANNOTATIONS_DIR}/gencode.v19.annotation.bed' > '${ANNOTATIONS_DIR}/gencode.v19.annotation.genes.bed' && \
            bedtools sort -i '${ANNOTATIONS_DIR}/gencode.v19.annotation.genes.bed' > '$SORTED_GENE_BED'
        "
        echo "    ✓ Annotations prepared"
    fi
else
    echo "[2/4] Skipping annotations (--skip-annotation)"
fi

# Step 3: Sample data (optional)
if [ -n "$SAMPLE_ID" ]; then
    FASTQ_1="data/inputs/fastq/${SAMPLE_ID}_1.fastq"
    if [ -f "$FASTQ_1" ]; then
        echo "[3/4] Sample $SAMPLE_ID already exists, skipping..."
    else
        echo "[3/4] Downloading sample: $SAMPLE_ID (this may take a while)..."
        mkdir -p data/inputs/fastq
        
        run_docker fasterq-dump "$SAMPLE_ID" --split-files -O data/inputs/fastq/
        echo "    ✓ Sample downloaded"
    fi
else
    echo "[3/4] Skipping sample download (use --sample SRR413984 to download)"
fi

# Step 4: BWA index
if [ "$SKIP_INDEX" = false ]; then
    if [ -f "$REFERENCE" ]; then
        if [ -f "${REFERENCE}.bwt.2bit.64" ]; then
            echo "[4/4] BWA index already exists, skipping..."
        else
            echo "[4/4] Creating BWA-MEM2 index (this takes ~1 hour for full genome)..."
            run_docker bwa-mem2 index "$REFERENCE"
            run_docker samtools faidx "$REFERENCE"
            echo "    ✓ BWA index created"
        fi
    else
        echo "[4/4] Skipping BWA index (no reference genome found)"
    fi
else
    echo "[4/4] Skipping BWA index (--skip-index)"
fi

echo ""
echo "========================================"
echo "  ✅ Data setup complete!"
echo "========================================"
echo ""
echo "Next steps:"
if [ -n "$SAMPLE_ID" ]; then
    echo "  ./docker/run_alignment.sh $SAMPLE_ID 8"
    echo "  ./docker/run_full_pipeline.sh data/intermediate/${SAMPLE_ID}.sorted.bam \\"
    echo "      $REFERENCE $SORTED_GENE_BED data/outputs/${SAMPLE_ID}"
fi
