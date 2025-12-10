#!/usr/bin/env python3
"""
get_seq.py - Extract and analyze sequences around eccDNA junctions.

Uses Smith-Waterman alignment to detect microhomology at eccDNA junction sites
by comparing clipped consensus sequences with reference genome sequences.
"""

import argparse
import sys
from pathlib import Path

import pysam

from . import sw


def analyze_junction(
    ref_genome: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    start_seq: str,
    end_seq: str,
    gap: int = -2,
    mismatch: int = -1,
    match: int = 1,
) -> None:
    """
    Analyze eccDNA junction for microhomology.

    Args:
        ref_genome: Opened pysam FastaFile.
        chrom: Chromosome name in the reference.
        start: Start position of the eccDNA candidate.
        end: End position of the eccDNA candidate.
        start_seq: Consensus sequence from soft clips at the start.
        end_seq: Consensus sequence from soft clips at the end.
        gap: Gap penalty for alignment.
        mismatch: Mismatch penalty for alignment.
        match: Match score for alignment.
    """
    # Fetch the reference sequence for this region
    seq = ref_genome.fetch(chrom, start, end)
    print(f"Reference sequence ({len(seq)} bp):")
    print(seq[:50] + "..." if len(seq) > 50 else seq)

    # Extract flanking sequences for microhomology detection
    flank_start = seq[:10]
    flank_end = seq[-10:]

    # Find microhomology between junction flanks
    scoring_matrix = sw.sw_fill_matrix(flank_start, flank_end, gap, mismatch, match)
    align_a, align_b, score = sw.sw_traceback(
        scoring_matrix, flank_start, flank_end, gap, mismatch, match
    )

    print("\nMicrohomology detection (flanks):")
    print(f"  5' flank:  {flank_start}")
    print(f"  3' flank:  {flank_end}")
    print(f"  Aligned A: {align_a}")
    print(f"  Aligned B: {align_b}")
    print(f"  Score:     {score}")

    # Build predicted circular junction sequence
    micro_h = align_a
    c_begin = micro_h + end_seq
    c_end = start_seq + micro_h

    # Align predicted junction with reference
    print(f"\nPredicted 5' junction: {c_begin}")
    junction_ref_start = seq[: len(c_begin)]
    h1 = sw.sw_fill_matrix(junction_ref_start, c_begin, gap, mismatch, match)
    a1, b1, s1 = sw.sw_traceback(h1, junction_ref_start, c_begin, gap, mismatch, match)
    print(f"  Ref:       {a1}")
    print(f"  Junction:  {b1}")
    print(f"  Score:     {s1}")

    print(f"\nPredicted 3' junction: {c_end}")
    junction_ref_end = seq[-len(c_end) :]
    h2 = sw.sw_fill_matrix(junction_ref_end, c_end, gap, mismatch, match)
    a2, b2, s2 = sw.sw_traceback(h2, junction_ref_end, c_end, gap, mismatch, match)
    print(f"  Ref:       {a2}")
    print(f"  Junction:  {b2}")
    print(f"  Score:     {s2}")


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Extract and analyze sequences around eccDNA junctions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--fasta", type=str, required=True,
        help="Path to the reference genome FASTA file",
    )
    parser.add_argument(
        "--chrom", type=str, required=True,
        help="Chromosome name in the reference",
    )
    parser.add_argument(
        "--start", type=int, required=True,
        help="Start position of the eccDNA candidate",
    )
    parser.add_argument(
        "--end", type=int, required=True,
        help="End position of the eccDNA candidate",
    )
    parser.add_argument(
        "--start_seq", type=str, required=True,
        help="Consensus sequence from soft clips at the start",
    )
    parser.add_argument(
        "--end_seq", type=str, required=True,
        help="Consensus sequence from soft clips at the end",
    )
    parser.add_argument(
        "--gap", type=int, default=-2,
        help="Gap penalty for alignment",
    )
    parser.add_argument(
        "--mismatch", type=int, default=-1,
        help="Mismatch penalty for alignment",
    )
    parser.add_argument(
        "--match", type=int, default=1,
        help="Match score for alignment",
    )
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        print(f"Error: FASTA file not found: {fasta_path}", file=sys.stderr)
        sys.exit(1)

    with pysam.FastaFile(args.fasta) as ref_genome:
        analyze_junction(
            ref_genome,
            args.chrom,
            args.start,
            args.end,
            args.start_seq,
            args.end_seq,
            args.gap,
            args.mismatch,
            args.match,
        )


if __name__ == "__main__":
    main()
