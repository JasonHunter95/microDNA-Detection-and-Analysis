#!/usr/bin/env python3
"""
microDNA.py - Detect microDNA candidates from soft-clipped reads in BAM files.

This script analyzes BAM files to identify potential microDNA (small extrachromosomal
circular DNA) by examining soft-clipped reads and building consensus sequences.
"""

import argparse
import sys
from collections import Counter
from typing import Optional

import pysam


def most_common_char(char_list: list[str]) -> Optional[str]:
    """
    Return the most common character in a list.

    Args:
        char_list: A list of single characters.

    Returns:
        The most frequently occurring character, or None if the list is empty.
    """
    if not char_list:
        return None
    count = Counter(char_list)
    return count.most_common(1)[0][0]


def get_clipped_seq(read: pysam.AlignedSegment) -> str:
    """
    Extract the soft-clipped sequence from a read.

    Properly handles all CIGAR operations that consume query bases:
    - 0 (M): Match/Mismatch
    - 1 (I): Insertion to reference
    - 4 (S): Soft clip
    - 7 (=): Sequence match
    - 8 (X): Sequence mismatch

    Args:
        read: A pysam AlignedSegment object.

    Returns:
        A string containing only the soft-clipped portions of the read sequence.
    """
    clipped_seq = ""
    if read.cigartuples is None:
        return clipped_seq

    # CIGAR operations that consume query sequence
    QUERY_CONSUMING_OPS = {0, 1, 4, 7, 8}  # M, I, S, =, X

    pos = 0
    for op, length in read.cigartuples:
        if op == 4:  # Soft clip - extract and advance
            clipped_seq += read.query_sequence[pos : pos + length]
            pos += length
        elif op in QUERY_CONSUMING_OPS:  # Other query-consuming ops - just advance
            pos += length
        # Operations 2 (D), 3 (N), 5 (H), 6 (P) don't consume query
    return clipped_seq


def count_matching_chars(char: str, char_list: list[str]) -> int:
    """Count occurrences of a character in a list."""
    return sum(1 for c in char_list if c == char)


def get_consensus(
    seq_list: list[str],
    min_num_chars: int = 15,
    min_char_density: float = 0.75,
) -> str:
    """
    Build a consensus sequence from a list of sequences.

    Args:
        seq_list: List of sequences to build consensus from.
        min_num_chars: Minimum number of characters at a position to include in consensus.
        min_char_density: Minimum density of the most common character to include.

    Returns:
        The consensus sequence string.
    """
    consensus: dict[int, list[str]] = {}
    for seq in seq_list:
        for i, char in enumerate(seq):
            if i not in consensus:
                consensus[i] = []
            consensus[i].append(char)

    result = ""
    for i in sorted(consensus.keys()):
        mcc = most_common_char(consensus[i])
        if mcc is None:
            continue
        mcc_count = count_matching_chars(mcc, consensus[i])
        density = mcc_count / len(consensus[i])
        if mcc_count >= min_num_chars and density >= min_char_density:
            result += mcc
    return result


def process_bam(
    bam_path: str,
    min_threshold: int = 50,
    max_threshold: int = 200,
    min_num_chars: int = 15,
    min_char_density: float = 0.75,
) -> None:
    """
    Process a BAM file to find microDNA candidates.

    Args:
        bam_path: Path to the input BAM file.
        min_threshold: Minimum number of reads in a window to report.
        max_threshold: Maximum number of reads in a window to report.
        min_num_chars: Minimum characters for consensus building.
        min_char_density: Minimum density for consensus building.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        start_window: list[tuple[int, str]] = []
        end_window: list[tuple[int, str]] = []

        for read in bam:
            if read.cigartuples is None:
                continue

            # Check for soft clip at the start of the read
            if read.cigartuples[0][0] == 4:
                if len(start_window) == 0:
                    start_window.append((read.pos, get_clipped_seq(read)))
                elif read.pos == start_window[0][0]:
                    start_window.append((read.pos, get_clipped_seq(read)))
                else:
                    if min_threshold < len(start_window) < max_threshold:
                        print(len(start_window), start_window[0], "+")
                        s_s = get_consensus(
                            [s[1][::-1] for s in start_window],
                            min_num_chars,
                            min_char_density,
                        )
                        print(s_s[::-1])
                    start_window = [(read.pos, get_clipped_seq(read))]

            # Check for soft clip at the end of the read
            elif read.cigartuples[-1][0] == 4:
                circle_end = read.pos + read.cigartuples[0][1]
                if len(end_window) == 0:
                    end_window.append((circle_end, get_clipped_seq(read)))
                elif circle_end == end_window[0][0]:
                    end_window.append((circle_end, get_clipped_seq(read)))
                else:
                    if min_threshold < len(end_window) < max_threshold:
                        print(len(end_window), end_window[0], "-")
                        s_e = get_consensus(
                            [s[1] for s in end_window],
                            min_num_chars,
                            min_char_density,
                        )
                        print(s_e)
                    end_window = [(circle_end, get_clipped_seq(read))]


def main() -> None:
    """Main entry point for the microDNA detection script."""
    parser = argparse.ArgumentParser(
        description="Detect microDNA candidates from soft-clipped reads in BAM files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the input sorted BAM file.",
    )
    parser.add_argument(
        "--min-threshold",
        type=int,
        default=50,
        help="Minimum number of reads in a window to report.",
    )
    parser.add_argument(
        "--max-threshold",
        type=int,
        default=200,
        help="Maximum number of reads in a window to report.",
    )
    parser.add_argument(
        "--min-num-chars",
        type=int,
        default=15,
        help="Minimum number of characters for consensus building.",
    )
    parser.add_argument(
        "--min-char-density",
        type=float,
        default=0.75,
        help="Minimum character density for consensus building.",
    )

    args = parser.parse_args()

    try:
        process_bam(
            args.input,
            args.min_threshold,
            args.max_threshold,
            args.min_num_chars,
            args.min_char_density,
        )
    except FileNotFoundError:
        print(f"Error: BAM file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing BAM file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
