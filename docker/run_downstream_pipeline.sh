#!/bin/bash
# run_downstream_pipeline.sh
# Runs downstream steps (clean_bed, bedtools closest, annotate) in Docker
#
# Usage:
#   ./docker/run_downstream_pipeline.sh <eccdna_bed> <annotation_bed> <output_prefix>
#
# Example:
#   ./docker/run_downstream_pipeline.sh \
#       data/intermediate/SRR413984.eccdna.bed \
#       data/inputs/references/annotations/gencode.v19.annotation.genes.sorted.bed \
#       data/outputs/SRR413984

set -e

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <eccdna_bed> <annotation_bed> <output_prefix>"
    echo ""
    echo "Arguments:"
    echo "  eccdna_bed      Path to raw eccDNA BED file from Circle-Map"
    echo "  annotation_bed  Path to sorted gene annotation BED file"
    echo "  output_prefix   Prefix for output files (e.g., data/outputs/SRR413984)"
    exit 1
fi

ECCDNA_BED="$1"
ANNOTATION_BED="$2"
OUTPUT_PREFIX="$3"

# Derived paths
CLEANED_BED="${OUTPUT_PREFIX}.eccdna.cleaned.bed"
SORTED_BED="${OUTPUT_PREFIX}.eccdna.sorted.bed"
CLOSEST_BED="${OUTPUT_PREFIX}.eccdna.closest_genes.bed"
ANNOTATED_TSV="${OUTPUT_PREFIX}.eccdna.annotated.tsv"

IMAGE_NAME="microdna-downstream:latest"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== microDNA Downstream Pipeline (Docker) ==="
echo "Input BED:     $ECCDNA_BED"
echo "Annotation:    $ANNOTATION_BED"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Build Docker image if needed
if ! docker image inspect "$IMAGE_NAME" &> /dev/null; then
    echo "Building Docker image..."
    docker build -t "$IMAGE_NAME" -f "$PROJECT_DIR/docker/Dockerfile.microdna" "$PROJECT_DIR/"
fi

# Helper function to run Python modules in Docker
run_docker() {
    docker run --rm -v "$PROJECT_DIR:/workspace" -w /workspace "$IMAGE_NAME" "$@"
}

# Helper function to run shell commands in Docker
run_docker_shell() {
    docker run --rm -v "$PROJECT_DIR:/workspace" -w /workspace --entrypoint bash "$IMAGE_NAME" -c "$1"
}

echo "[1/4] Cleaning BED format..."
run_docker microdna.clean_bed -i "$ECCDNA_BED" -o "$CLEANED_BED"

echo "[2/4] Sorting eccDNA BED..."
run_docker_shell "bedtools sort -i $CLEANED_BED > $SORTED_BED"

echo "[3/4] Finding closest genes..."
run_docker_shell "bedtools closest -a $SORTED_BED -b $ANNOTATION_BED -d > $CLOSEST_BED"

echo "[4/4] Generating annotated output..."
run_docker microdna.annotate_eccdna_closest -i "$CLOSEST_BED" -o "$ANNOTATED_TSV"

echo ""
echo "✅ Downstream pipeline complete!"
echo "Output: $ANNOTATED_TSV"
