#!/usr/bin/env python3
"""
get_soft_clips.py - Count soft-clipped reads in a BAM file.

Analyzes a BAM file to count reads with soft clips at the start and/or end,
which are potential indicators of circular DNA junction reads.
"""

import argparse
import sys
from pathlib import Path

import pysam


def count_soft_clips(bam_path: str) -> dict[str, int]:
    """
    Count soft-clipped reads in a BAM file.

    Args:
        bam_path: Path to the input BAM file.

    Returns:
        Dictionary with counts for total_reads, start_softclip, end_softclip.
    """
    start_softclip = 0
    end_softclip = 0
    total_reads = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            total_reads += 1
            if read.cigartuples is None:
                continue
            # Soft-clip at start (CIGAR op 4)
            if read.cigartuples[0][0] == 4:
                start_softclip += 1
            # Soft-clip at end (CIGAR op 4)
            if read.cigartuples[-1][0] == 4:
                end_softclip += 1

    return {
        "total_reads": total_reads,
        "start_softclip": start_softclip,
        "end_softclip": end_softclip,
    }


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Count soft-clipped reads in a BAM file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to the input BAM file",
    )
    args = parser.parse_args()

    bam_path = Path(args.input)
    if not bam_path.exists():
        print(f"Error: BAM file not found: {bam_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Analyzing: {args.input}")
    counts = count_soft_clips(args.input)

    print("\nResults:")
    print(f"  Total reads:             {counts['total_reads']:,}")
    print(f"  Soft-clipped at start:   {counts['start_softclip']:,}")
    print(f"  Soft-clipped at end:     {counts['end_softclip']:,}")
    print(f"  Total soft-clipped:      {counts['start_softclip'] + counts['end_softclip']:,}")


if __name__ == "__main__":
    main()
