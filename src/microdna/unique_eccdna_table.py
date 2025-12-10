#!/usr/bin/env python3
"""
unique_eccdna_table.py - Generate a unique eccDNA summary table.

This script aggregates eccDNA annotations to produce a table with one row per
unique eccDNA, summarizing overlapping gene information.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def aggregate_eccdna(input_path: str, output_path: str) -> None:
    """
    Aggregate eccDNA annotations into a unique eccDNA summary table.

    Args:
        input_path: Path to input TSV with eccDNA-gene intersection annotations.
        output_path: Path to output unique eccDNA summary TSV.
    """
    df_eccdna = pd.read_csv(input_path, sep="\t", header=0)

    eccdna_unique = df_eccdna.groupby("ecc_id").agg(
        ecc_chr=("ecc_chr", "first"),
        ecc_start=("ecc_start", "first"),
        ecc_end=("ecc_end", "first"),
        ecc_length=("ecc_length", "first"),
        ecc_score=("ecc_score", "first"),
        num_overlapping_genes=("gene_id", "nunique"),
        overlapping_gene_ids=("gene_id", lambda x: list(x.dropna().unique())),
        overlapping_gene_names=("gene_name", lambda x: list(x.dropna().unique())),
        overlapping_gene_types=("gene_type", lambda x: list(x.dropna().unique())),
    ).reset_index()

    eccdna_unique.to_csv(output_path, sep="\t", index=False)

    print(f"✅ Unique eccDNA table ({len(eccdna_unique)} rows) saved to: {output_path}")


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate a unique eccDNA summary table from intersection annotations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input TSV with eccDNA-gene intersection annotations",
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output unique eccDNA summary TSV",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    aggregate_eccdna(args.input, args.output)


if __name__ == "__main__":
    main()