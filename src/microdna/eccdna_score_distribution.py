#!/usr/bin/env python3
"""
eccdna_score_distribution.py - Plot eccDNA confidence score distribution.

Generates a histogram with KDE overlay showing the distribution of eccDNA
confidence scores, including statistical markers (mean, median, quartiles).
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from .plotting import plot_distribution


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Plot eccDNA score distribution from a TSV file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to input eccDNA TSV file with ecc_score column",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Path to output PNG file",
    )
    parser.add_argument(
        "--show", action="store_true",
        help="Display the plot interactively",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    ecc_data = pd.read_csv(args.input, sep="\t", low_memory=False)

    if "ecc_score" not in ecc_data.columns:
        print("Error: Required column 'ecc_score' missing from the data.", file=sys.stderr)
        sys.exit(1)

    plot_distribution(
        ecc_data,
        "ecc_score",
        "Distribution of eccDNA Confidence Scores (Log Scale)",
        "Score [Log Scale]",
    )

    # Ensure output directory exists
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.output, dpi=150, bbox_inches="tight")
    print(f"✅ Saved plot to: {args.output}")

    if args.show:
        plt.show()

    plt.close()


if __name__ == "__main__":
    main()