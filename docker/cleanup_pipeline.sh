#!/bin/bash
# cleanup_pipeline.sh - Remove intermediate pipeline files to reclaim disk space
#
# This script is designed to be run after the pipeline has completed successfully
# and you have verified the final output files. It removes large intermediate files
# that are not needed for downstream analysis.
#
# IMPORTANT: Always verify your final outputs before running cleanup!
#
# Usage:
#   ./docker/cleanup_pipeline.sh <sample_id> [options]
#
# Options:
#   --dry-run       Show what would be deleted without deleting (default)
#   --execute       Actually delete the files
#   --keep-bam      Keep the aligned BAM file (useful for future reanalysis)
#   --keep-fastq    Keep the raw FASTQ files
#   --all           Remove everything including FASTQ and aligned BAM
#
# Examples:
#   # Preview what would be deleted (safe)
#   ./docker/cleanup_pipeline.sh SRR413984
#
#   # Delete intermediate files but keep aligned BAM and FASTQ
#   ./docker/cleanup_pipeline.sh SRR413984 --execute
#
#   # Delete everything including FASTQ and aligned BAM
#   ./docker/cleanup_pipeline.sh SRR413984 --execute --all

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

usage() {
    echo "Usage: $0 <sample_id> [options]"
    echo ""
    echo "Options:"
    echo "  --dry-run       Show what would be deleted without deleting (default)"
    echo "  --execute       Actually delete the files"
    echo "  --keep-bam      Keep the aligned BAM file (default)"
    echo "  --keep-fastq    Keep the raw FASTQ files (default)"
    echo "  --all           Remove everything including FASTQ and aligned BAM"
    echo ""
    echo "Example:"
    echo "  $0 SRR413984 --execute"
    exit 1
}

if [ "$#" -lt 1 ]; then
    usage
fi

SAMPLE_ID="$1"
shift

# Default options
DRY_RUN=true
KEEP_BAM=true
KEEP_FASTQ=true

# Parse options
while [ "$#" -gt 0 ]; do
    case "$1" in
        --dry-run)
            DRY_RUN=true
            ;;
        --execute)
            DRY_RUN=false
            ;;
        --keep-bam)
            KEEP_BAM=true
            ;;
        --keep-fastq)
            KEEP_FASTQ=true
            ;;
        --all)
            KEEP_BAM=false
            KEEP_FASTQ=false
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
    shift
done

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

# Define file categories with their paths
declare -a UPSTREAM_INTERMEDIATES=(
    "data/intermediate/${SAMPLE_ID}.querysorted.bam"
    "data/intermediate/${SAMPLE_ID}.candidates.bam"
    "data/intermediate/${SAMPLE_ID}.candidates.sorted.bam"
    "data/intermediate/${SAMPLE_ID}.candidates.sorted.bam.bai"
)

declare -a DOWNSTREAM_INTERMEDIATES=(
    "data/outputs/${SAMPLE_ID}.eccdna.bed"
    "data/outputs/${SAMPLE_ID}.eccdna.cleaned.bed"
    "data/outputs/${SAMPLE_ID}.eccdna.sorted.bed"
    "data/outputs/${SAMPLE_ID}.eccdna.closest_genes.bed"
)

declare -a ALIGNMENT_FILES=(
    "data/intermediate/${SAMPLE_ID}.sorted.bam"
    "data/intermediate/${SAMPLE_ID}.sorted.bam.bai"
)

declare -a FASTQ_FILES=(
    "data/inputs/fastq/${SAMPLE_ID}_1.fastq"
    "data/inputs/fastq/${SAMPLE_ID}_2.fastq"
    "data/inputs/fastq/${SAMPLE_ID}_1.fastq.gz"
    "data/inputs/fastq/${SAMPLE_ID}_2.fastq.gz"
)

# Calculate sizes and track what to delete
calculate_size() {
    local file="$1"
    if [ -f "$file" ]; then
        # macOS compatible size calculation
        stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo 0
    else
        echo 0
    fi
}

human_readable_size() {
    local bytes="$1"
    if [ "$bytes" -ge 1073741824 ]; then
        echo "$(echo "scale=2; $bytes / 1073741824" | bc)GB"
    elif [ "$bytes" -ge 1048576 ]; then
        echo "$(echo "scale=2; $bytes / 1048576" | bc)MB"
    elif [ "$bytes" -ge 1024 ]; then
        echo "$(echo "scale=2; $bytes / 1024" | bc)KB"
    else
        echo "${bytes}B"
    fi
}

echo ""
echo "=========================================="
echo "  Pipeline Cleanup - Sample: $SAMPLE_ID"
echo "=========================================="

if [ "$DRY_RUN" = true ]; then
    echo -e "${YELLOW}MODE: DRY RUN (no files will be deleted)${NC}"
else
    echo -e "${RED}MODE: EXECUTE (files will be permanently deleted!)${NC}"
fi
echo ""

# Check final output exists before allowing cleanup
FINAL_OUTPUT="data/outputs/${SAMPLE_ID}.eccdna.annotated.tsv"
if [ ! -f "$FINAL_OUTPUT" ]; then
    echo -e "${RED}Warning: Final output not found: $FINAL_OUTPUT${NC}"
    echo "Run the full pipeline before cleanup, or verify sample ID."
    echo ""
fi

TOTAL_SIZE=0
FILES_TO_DELETE=()

echo "--- Upstream Intermediates ---"
for file in "${UPSTREAM_INTERMEDIATES[@]}"; do
    if [ -f "$file" ]; then
        size=$(calculate_size "$file")
        TOTAL_SIZE=$((TOTAL_SIZE + size))
        FILES_TO_DELETE+=("$file")
        echo -e "  ${RED}DELETE${NC}: $file ($(human_readable_size $size))"
    fi
done

echo ""
echo "--- Downstream Intermediates ---"
for file in "${DOWNSTREAM_INTERMEDIATES[@]}"; do
    if [ -f "$file" ]; then
        size=$(calculate_size "$file")
        TOTAL_SIZE=$((TOTAL_SIZE + size))
        FILES_TO_DELETE+=("$file")
        echo -e "  ${RED}DELETE${NC}: $file ($(human_readable_size $size))"
    fi
done

echo ""
echo "--- Aligned BAM Files ---"
for file in "${ALIGNMENT_FILES[@]}"; do
    if [ -f "$file" ]; then
        size=$(calculate_size "$file")
        if [ "$KEEP_BAM" = true ]; then
            echo -e "  ${GREEN}KEEP${NC}:   $file ($(human_readable_size $size))"
        else
            TOTAL_SIZE=$((TOTAL_SIZE + size))
            FILES_TO_DELETE+=("$file")
            echo -e "  ${RED}DELETE${NC}: $file ($(human_readable_size $size))"
        fi
    fi
done

echo ""
echo "--- Raw FASTQ Files ---"
for file in "${FASTQ_FILES[@]}"; do
    if [ -f "$file" ]; then
        size=$(calculate_size "$file")
        if [ "$KEEP_FASTQ" = true ]; then
            echo -e "  ${GREEN}KEEP${NC}:   $file ($(human_readable_size $size))"
        else
            TOTAL_SIZE=$((TOTAL_SIZE + size))
            FILES_TO_DELETE+=("$file")
            echo -e "  ${RED}DELETE${NC}: $file ($(human_readable_size $size))"
        fi
    fi
done

echo ""
echo "--- Final Outputs (always kept) ---"
for file in "data/outputs/${SAMPLE_ID}.eccdna.annotated.tsv"; do
    if [ -f "$file" ]; then
        size=$(calculate_size "$file")
        echo -e "  ${BLUE}FINAL${NC}:  $file ($(human_readable_size $size))"
    fi
done

echo ""
echo "=========================================="
echo -e "Total to reclaim: ${GREEN}$(human_readable_size $TOTAL_SIZE)${NC}"
echo "Files to delete:  ${#FILES_TO_DELETE[@]}"
echo "=========================================="
echo ""

if [ ${#FILES_TO_DELETE[@]} -eq 0 ]; then
    echo "No files to delete."
    exit 0
fi

if [ "$DRY_RUN" = true ]; then
    echo -e "${YELLOW}This was a dry run. To actually delete files, run:${NC}"
    echo "  $0 $SAMPLE_ID --execute"
    if [ "$KEEP_BAM" = true ] || [ "$KEEP_FASTQ" = true ]; then
        echo ""
        echo "To also remove FASTQ and aligned BAM files:"
        echo "  $0 $SAMPLE_ID --execute --all"
    fi
else
    echo -e "${RED}Deleting files...${NC}"
    for file in "${FILES_TO_DELETE[@]}"; do
        rm -f "$file"
        echo "  Deleted: $file"
    done
    echo ""
    echo -e "${GREEN}✅ Cleanup complete! Reclaimed $(human_readable_size $TOTAL_SIZE)${NC}"
fi
