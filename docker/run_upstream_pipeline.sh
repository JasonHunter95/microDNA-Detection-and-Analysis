#!/bin/bash
# run_upstream_pipeline.sh
# Runs Circle-Map upstream steps (ReadExtractor + Realign) in Docker
#
# Usage:
#   ./docker/run_upstream_pipeline.sh <sorted_bam> <reference_fasta> <output_prefix>
#
# Example:
#   ./docker/run_upstream_pipeline.sh \
#       data/intermediate/SRR413984.sorted.bam \
#       data/inputs/references/genome/GCF_000001405.13_GRCh37_genomic.fna \
#       data/intermediate/SRR413984

set -e

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <sorted_bam> <reference_fasta> <output_prefix>"
    echo ""
    echo "Arguments:"
    echo "  sorted_bam      Path to position-sorted BAM file (relative to project root)"
    echo "  reference_fasta Path to reference FASTA file (relative to project root)"
    echo "  output_prefix   Prefix for output files (e.g., data/intermediate/SRR413984)"
    exit 1
fi

SORTED_BAM="$1"
REFERENCE="$2"
OUTPUT_PREFIX="$3"

# Derive paths
QUERYSORTED_BAM="${OUTPUT_PREFIX}.querysorted.bam"
CANDIDATES_BAM="${OUTPUT_PREFIX}.candidates.bam"
CANDIDATES_SORTED_BAM="${OUTPUT_PREFIX}.candidates.sorted.bam"
ECCDNA_BED="${OUTPUT_PREFIX}.eccdna.bed"

IMAGE_NAME="microdna-circlemap:latest"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== microDNA Upstream Pipeline (Docker) ==="
echo "Input BAM:     $SORTED_BAM"
echo "Reference:     $REFERENCE"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Build Docker image if needed
if ! docker image inspect "$IMAGE_NAME" &> /dev/null; then
    echo "Building Docker image..."
    docker build -t "$IMAGE_NAME" -f "$PROJECT_DIR/docker/Dockerfile.circlemap" "$PROJECT_DIR/docker/"
fi

# Helper function to run commands in Docker
# Mounts project as /workspace and sets working directory
run_docker() {
    docker run --rm -v "$PROJECT_DIR:/workspace" -w /workspace "$IMAGE_NAME" "$@"
}

echo "[1/5] Query-sorting BAM..."
run_docker samtools sort -n -o "$QUERYSORTED_BAM" "$SORTED_BAM"

echo "[2/5] Extracting candidate circular reads..."
run_docker Circle-Map ReadExtractor \
    -i "$QUERYSORTED_BAM" \
    -o "$CANDIDATES_BAM"

echo "[3/5] Sorting candidate BAM..."
run_docker samtools sort -o "$CANDIDATES_SORTED_BAM" "$CANDIDATES_BAM"

echo "[4/5] Indexing candidate BAM..."
run_docker samtools index "$CANDIDATES_SORTED_BAM"

echo "[5/5] Realigning to detect eccDNA..."
run_docker Circle-Map Realign \
    -i "$CANDIDATES_SORTED_BAM" \
    -qbam "$QUERYSORTED_BAM" \
    -sbam "$SORTED_BAM" \
    -fasta "$REFERENCE" \
    -o "$ECCDNA_BED" \
    --split 2

echo ""
echo "✅ Upstream pipeline complete!"
echo "Output: $ECCDNA_BED"
echo ""
echo "Next steps (run locally):"
echo "  python -m microdna.clean_bed -i $ECCDNA_BED -o ${OUTPUT_PREFIX}.eccdna.cleaned.bed"
