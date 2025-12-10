#!/usr/bin/env python3
"""
sw.py - Smith-Waterman local sequence alignment algorithm.

Implements the Smith-Waterman algorithm for local sequence alignment,
commonly used in bioinformatics to identify regions of similarity between
DNA or protein sequences. This is particularly useful for detecting
microhomology at eccDNA junction sites.

The algorithm uses dynamic programming to find the optimal local alignment,
allowing for gaps and mismatches with configurable scoring parameters.
"""

import argparse


def sw_fill_matrix(
    seq_a: str,
    seq_b: str,
    gap: int,
    mismatch: int,
    match: int,
) -> list[list[int]]:
    """
    Fill the Smith-Waterman scoring matrix.

    Args:
        seq_a: First sequence to align.
        seq_b: Second sequence to align.
        gap: Penalty for introducing a gap (typically negative).
        mismatch: Penalty for a mismatch (typically negative).
        match: Score for a match (typically positive).

    Returns:
        2D scoring matrix where H[i][j] represents the best local alignment
        score ending at position i in seq_a and position j in seq_b.
    """
    rows = len(seq_a) + 1
    cols = len(seq_b) + 1
    scoring_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    for i in range(1, rows):
        for j in range(1, cols):
            # Score for match/mismatch at current position
            diag_score = scoring_matrix[i - 1][j - 1] + (
                match if seq_a[i - 1] == seq_b[j - 1] else mismatch
            )
            # Score for gap in seq_b (vertical move)
            up_score = scoring_matrix[i - 1][j] + gap
            # Score for gap in seq_a (horizontal move)
            left_score = scoring_matrix[i][j - 1] + gap

            # Local alignment: reset to 0 if all paths are negative
            scoring_matrix[i][j] = max(diag_score, up_score, left_score, 0)

    return scoring_matrix


def sw_traceback(
    scoring_matrix: list[list[int]],
    seq_a: str,
    seq_b: str,
    gap: int,
    mismatch: int,
    match: int,
) -> tuple[str, str, int]:
    """
    Perform traceback to find the optimal local alignment.

    Args:
        scoring_matrix: Filled Smith-Waterman scoring matrix.
        seq_a: First sequence.
        seq_b: Second sequence.
        gap: Gap penalty.
        mismatch: Mismatch penalty.
        match: Match score.

    Returns:
        Tuple of (aligned_seq_a, aligned_seq_b, alignment_score).
    """
    # Find the cell with the maximum score
    max_score = 0
    max_i, max_j = 0, 0
    for i in range(len(seq_a) + 1):
        for j in range(len(seq_b) + 1):
            if scoring_matrix[i][j] > max_score:
                max_score = scoring_matrix[i][j]
                max_i, max_j = i, j

    # Traceback from maximum score cell
    i, j = max_i, max_j
    align_a: list[str] = []
    align_b: list[str] = []

    while scoring_matrix[i][j] > 0:
        current_score = scoring_matrix[i][j]
        diag_score = scoring_matrix[i - 1][j - 1] + (
            match if seq_a[i - 1] == seq_b[j - 1] else mismatch
        )

        if current_score == diag_score:
            # Diagonal move (match/mismatch)
            align_a.append(seq_a[i - 1])
            align_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif current_score == scoring_matrix[i - 1][j] + gap:
            # Vertical move (gap in seq_b)
            align_a.append(seq_a[i - 1])
            align_b.append("-")
            i -= 1
        elif current_score == scoring_matrix[i][j - 1] + gap:
            # Horizontal move (gap in seq_a)
            align_a.append("-")
            align_b.append(seq_b[j - 1])
            j -= 1
        else:
            break

    # Reverse alignments (traceback goes backwards)
    return "".join(align_a[::-1]), "".join(align_b[::-1]), max_score


def main() -> None:
    """Main entry point for Smith-Waterman alignment."""
    parser = argparse.ArgumentParser(
        description="Perform Smith-Waterman local sequence alignment.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--A", type=str, required=True,
        help="First sequence to align",
    )
    parser.add_argument(
        "--B", type=str, required=True,
        help="Second sequence to align",
    )
    parser.add_argument(
        "--gap", type=int, default=-2,
        help="Gap penalty (typically negative)",
    )
    parser.add_argument(
        "--miss", type=int, default=-1,
        help="Mismatch penalty (typically negative)",
    )
    parser.add_argument(
        "--match", type=int, default=1,
        help="Match score (typically positive)",
    )
    args = parser.parse_args()

    scoring_matrix = sw_fill_matrix(args.A, args.B, args.gap, args.miss, args.match)
    align_a, align_b, score = sw_traceback(
        scoring_matrix, args.A, args.B, args.gap, args.miss, args.match
    )

    print(f"Sequence A: {align_a}")
    print(f"Sequence B: {align_b}")
    print(f"Score: {score}")


if __name__ == "__main__":
    main()
