#!/bin/bash
# run_full_pipeline.sh
# Runs the complete microDNA detection pipeline in Docker (upstream + downstream)
#
# Usage:
#   ./docker/run_full_pipeline.sh <sorted_bam> <reference_fasta> <annotation_bed> <output_prefix>
#
# Example:
#   ./docker/run_full_pipeline.sh \
#       data/intermediate/SRR413984.sorted.bam \
#       data/inputs/references/genome/GCF_000001405.13_GRCh37_genomic.fna \
#       data/inputs/references/annotations/gencode.v19.annotation.genes.sorted.bed \
#       data/outputs/SRR413984

set -e

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <sorted_bam> <reference_fasta> <annotation_bed> <output_prefix>"
    echo ""
    echo "Arguments:"
    echo "  sorted_bam      Path to position-sorted BAM file"
    echo "  reference_fasta Path to reference FASTA file"
    echo "  annotation_bed  Path to sorted gene annotation BED file"
    echo "  output_prefix   Prefix for output files (e.g., data/outputs/SRR413984)"
    exit 1
fi

SORTED_BAM="$1"
REFERENCE="$2"
ANNOTATION_BED="$3"
OUTPUT_PREFIX="$4"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "========================================"
echo "  microDNA Full Pipeline (Docker)"
echo "========================================"
echo ""

echo "=== Phase 1: Upstream Pipeline ==="
"$SCRIPT_DIR/run_upstream_pipeline.sh" "$SORTED_BAM" "$REFERENCE" "$OUTPUT_PREFIX"

echo ""
echo "=== Phase 2: Downstream Pipeline ==="
"$SCRIPT_DIR/run_downstream_pipeline.sh" \
    "${OUTPUT_PREFIX}.eccdna.bed" \
    "$ANNOTATION_BED" \
    "$OUTPUT_PREFIX"

echo ""
echo "========================================"
echo "  ✅ Full pipeline complete!"
echo "========================================"
echo "Final output: ${OUTPUT_PREFIX}.eccdna.annotated.tsv"
echo ""
echo "To reclaim disk space after verifying results:"
echo "  ./docker/cleanup_pipeline.sh $(basename "$OUTPUT_PREFIX") --execute"
