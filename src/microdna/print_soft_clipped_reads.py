#!/usr/bin/env python3
"""
print_soft_clipped_reads.py - Display details of soft-clipped reads.

Prints detailed information about reads containing soft clips, useful for
debugging and understanding circular DNA junction signatures.
"""

import argparse
import sys
from pathlib import Path

import pysam


def print_soft_clipped_reads(bam_path: str, num_reads: int = 10) -> None:
    """
    Print details of soft-clipped reads from a BAM file.

    Args:
        bam_path: Path to the input BAM file.
        num_reads: Maximum number of reads to display.
    """
    count = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.cigartuples is None:
                continue

            # Check for soft clip at start or end
            has_start_clip = read.cigartuples[0][0] == 4
            has_end_clip = read.cigartuples[-1][0] == 4

            if has_start_clip or has_end_clip:
                clip_location = []
                if has_start_clip:
                    clip_location.append("start")
                if has_end_clip:
                    clip_location.append("end")

                print(f"--- Read {count + 1} ---")
                print(f"  Name:     {read.query_name}")
                print(f"  Position: {read.reference_name}:{read.reference_start}")
                print(f"  CIGAR:    {read.cigarstring}")
                print(f"  MAPQ:     {read.mapping_quality}")
                print(f"  Flags:    {read.flag}")
                print(f"  Clip at:  {', '.join(clip_location)}")
                print(f"  Sequence: {read.query_sequence[:50]}...")
                print()

                count += 1
                if count >= num_reads:
                    break

    print(f"Displayed {count} soft-clipped reads.")


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Display details of soft-clipped reads from a BAM file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to the input BAM file",
    )
    parser.add_argument(
        "-n", "--num-reads", type=int, default=10,
        help="Maximum number of reads to display",
    )
    args = parser.parse_args()

    bam_path = Path(args.input)
    if not bam_path.exists():
        print(f"Error: BAM file not found: {bam_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Scanning: {args.input}\n")
    print_soft_clipped_reads(args.input, args.num_reads)


if __name__ == "__main__":
    main()