#!/usr/bin/env python3
"""
cleaner_tsv.py - Clean and filter eccDNA annotation data.

This script validates, cleans, and filters eccDNA intersection data,
applying quality thresholds and computing overlap percentages.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


DTYPE_SPEC = {
    "ecc_start": "int64",
    "ecc_end": "int64",
    "ecc_length": "int64",
    "ecc_score": "float64",
    "gene_start": "int64",
    "gene_end": "int64",
    "overlap_bp": "int64",
    "ecc_chr": "category",
    "gene_chr": "category",
    "gene_strand": "category",
    "gene_type": "category",
}


def load_data(input_path: str, chunk_size: int = 100000) -> pd.DataFrame:
    """
    Load TSV data with chunking for memory efficiency.

    Args:
        input_path: Path to input TSV file.
        chunk_size: Number of rows per chunk.

    Returns:
        Loaded DataFrame.
    """
    chunks = []
    try:
        for chunk in pd.read_csv(
            input_path, sep="\t", dtype=DTYPE_SPEC,
            chunksize=chunk_size, low_memory=False
        ):
            chunks.append(chunk)
        return pd.concat(chunks, ignore_index=True)
    except Exception:
        # Fallback: load without specified dtypes
        print("Loading without initial dtypes...")
        chunks = []
        for chunk in pd.read_csv(
            input_path, sep="\t", chunksize=chunk_size, low_memory=False
        ):
            chunks.append(chunk)
        return pd.concat(chunks, ignore_index=True)


def validate_and_clean(df: pd.DataFrame) -> pd.DataFrame:
    """
    Validate coordinates and clean the DataFrame.

    Args:
        df: Input DataFrame.

    Returns:
        Cleaned DataFrame with validation report printed.
    """
    # Convert numeric columns
    numeric_cols = [
        "ecc_start", "ecc_end", "ecc_length", "ecc_score",
        "gene_start", "gene_end", "overlap_bp"
    ]
    for col in numeric_cols:
        if col in df.columns and df[col].dtype == "object":
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Drop rows missing critical columns
    critical_cols = ["ecc_chr", "ecc_start", "ecc_end", "ecc_id", "gene_id", "overlap_bp"]
    existing_critical = [c for c in critical_cols if c in df.columns]
    df.dropna(subset=existing_critical, inplace=True)

    # Remove duplicates
    initial_rows = len(df)
    df.drop_duplicates(inplace=True)
    print(f"Removed {initial_rows - len(df)} duplicate rows.")

    # Validation report
    print("\nValidation Report:")
    invalid_ecc = df[df["ecc_start"] > df["ecc_end"]].shape[0]
    invalid_gene = df[df["gene_start"] > df["gene_end"]].shape[0]
    print(f"  Invalid ecc coords (start > end): {invalid_ecc}")
    print(f"  Invalid gene coords (start > end): {invalid_gene}")

    return df


def filter_and_annotate(
    df: pd.DataFrame,
    min_score: float = 10.0,
    min_overlap_bp: int = 50,
    min_overlap_pct: float = 10.0,
) -> pd.DataFrame:
    """
    Filter by score and compute overlap percentages.

    Args:
        df: Cleaned DataFrame.
        min_score: Minimum eccDNA score threshold.
        min_overlap_bp: Minimum overlap in base pairs.
        min_overlap_pct: Minimum overlap percentage of eccDNA.

    Returns:
        Filtered DataFrame with overlap annotations.
    """
    # Filter by score
    df_filtered = df[df["ecc_score"] >= min_score].copy()
    print(f"Rows after score filter (>= {min_score}): {len(df_filtered)}")

    # Compute overlap percentages
    df_filtered["gene_length"] = df_filtered["gene_end"] - df_filtered["gene_start"] + 1
    df_filtered["perc_overlap_ecc"] = df_filtered.apply(
        lambda row: (row["overlap_bp"] / row["ecc_length"] * 100)
        if row["ecc_length"] > 0 else 0,
        axis=1,
    )
    df_filtered["perc_overlap_gene"] = df_filtered.apply(
        lambda row: (row["overlap_bp"] / row["gene_length"] * 100)
        if pd.notna(row["gene_length"]) and row["gene_length"] > 0 else 0,
        axis=1,
    )

    # Filter for meaningful overlap
    df_meaningful = df_filtered[
        (df_filtered["overlap_bp"] >= min_overlap_bp) &
        (df_filtered["perc_overlap_ecc"] >= min_overlap_pct)
    ].copy()
    print(f"Rows after overlap filter: {len(df_meaningful)}")

    return df_meaningful


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Clean and filter eccDNA annotation data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to input eccDNA intersection TSV",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Path to output cleaned TSV",
    )
    parser.add_argument(
        "--min-score", type=float, default=10.0,
        help="Minimum eccDNA confidence score",
    )
    parser.add_argument(
        "--min-overlap-bp", type=int, default=50,
        help="Minimum overlap in base pairs",
    )
    parser.add_argument(
        "--min-overlap-pct", type=float, default=10.0,
        help="Minimum overlap percentage of eccDNA",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loading: {args.input}")
    df = load_data(args.input)
    print(f"Initial shape: {df.shape}")

    df = validate_and_clean(df)
    df = filter_and_annotate(
        df,
        min_score=args.min_score,
        min_overlap_bp=args.min_overlap_bp,
        min_overlap_pct=args.min_overlap_pct,
    )

    df.to_csv(args.output, sep="\t", index=False)
    print(f"\n✅ Cleaned data saved to: {args.output}")


if __name__ == "__main__":
    main()