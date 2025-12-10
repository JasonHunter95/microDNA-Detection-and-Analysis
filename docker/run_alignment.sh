#!/bin/bash
# run_alignment.sh - Align FASTQ reads to reference genome (Docker)
#
# Usage:
#   ./docker/run_alignment.sh <sample_id> [threads]
#
# Example:
#   ./docker/run_alignment.sh SRR413984 8

set -e

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <sample_id> [threads]"
    echo ""
    echo "Arguments:"
    echo "  sample_id   Sample ID (e.g., SRR413984)"
    echo "  threads     Number of threads (default: 4)"
    exit 1
fi

SAMPLE_ID="$1"
THREADS="${2:-4}"

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

IMAGE_NAME="microdna-circlemap:latest"
REFERENCE="data/inputs/references/genome/GCF_000001405.13_GRCh37_genomic.fna"
FASTQ_1="data/inputs/fastq/${SAMPLE_ID}_1.fastq"
FASTQ_2="data/inputs/fastq/${SAMPLE_ID}_2.fastq"
OUTPUT_BAM="data/intermediate/${SAMPLE_ID}.sorted.bam"

echo "=== Alignment Pipeline (Docker) ==="
echo "Sample:    $SAMPLE_ID"
echo "Threads:   $THREADS"
echo "Reference: $REFERENCE"
echo ""

# Validate inputs
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference genome not found: $REFERENCE"
    echo "Run: ./docker/setup_data.sh"
    exit 1
fi

if [ ! -f "$FASTQ_1" ] || [ ! -f "$FASTQ_2" ]; then
    echo "Error: FASTQ files not found for $SAMPLE_ID"
    echo "Run: ./docker/setup_data.sh --sample $SAMPLE_ID"
    exit 1
fi

if [ ! -f "${REFERENCE}.bwt.2bit.64" ]; then
    echo "Error: BWA index not found. Run: ./docker/setup_data.sh"
    exit 1
fi

# Build image if needed
if ! docker image inspect "$IMAGE_NAME" &> /dev/null; then
    echo "Building Docker image..."
    docker build -t "$IMAGE_NAME" -f "$PROJECT_DIR/docker/Dockerfile.circlemap" "$PROJECT_DIR/docker/"
fi

mkdir -p data/intermediate

# Run alignment in Docker
echo "[1/2] Aligning reads with BWA-MEM2..."
docker run --rm -v "$PROJECT_DIR:/workspace" -w /workspace "$IMAGE_NAME" \
    bash -c "bwa-mem2 mem -t $THREADS '$REFERENCE' '$FASTQ_1' '$FASTQ_2' | \
             samtools sort -@ $THREADS -m 1G -o '$OUTPUT_BAM'"

echo "[2/2] Indexing BAM..."
docker run --rm -v "$PROJECT_DIR:/workspace" -w /workspace "$IMAGE_NAME" \
    samtools index "$OUTPUT_BAM"

echo ""
echo "✅ Alignment complete!"
echo "Output: $OUTPUT_BAM"
echo ""
echo "Next step:"
echo "  ./docker/run_full_pipeline.sh $OUTPUT_BAM \\"
echo "      $REFERENCE data/inputs/references/annotations/gencode.v19.annotation.genes.sorted.bed \\"
echo "      data/outputs/$SAMPLE_ID"
