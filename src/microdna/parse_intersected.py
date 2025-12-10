#!/usr/bin/env python3
"""
parse_intersected.py - Parse bedtools intersect output for eccDNA-gene overlaps.

This script parses the output of `bedtools intersect -wa -wb -wo` to extract
eccDNA-gene overlap annotations, including gene metadata from GTF attributes.
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


FINAL_COLUMNS = [
    "ecc_chr", "ecc_start", "ecc_end", "ecc_id", "ecc_length",
    "ecc_score", "ecc_strand", "coverage", "split_reads", "discordant_mates",
    "gene_id", "gene_name", "gene_type",
    "gene_chr", "gene_start", "gene_end", "gene_strand",
    "overlap_bp",
]


def parse_attributes(attr_string: str) -> dict[str, str]:
    """
    Parse GTF-style attribute string to extract gene identifiers.

    Args:
        attr_string: GTF attribute string (key "value"; format).

    Returns:
        Dictionary with gene_id, gene_name, and gene_type.
    """
    if not isinstance(attr_string, str):
        return {"gene_id": ".", "gene_name": ".", "gene_type": "."}

    attributes = {}
    for item in attr_string.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(" ", 1)
        if len(parts) == 2:
            key = parts[0].strip()
            value = parts[1].strip()
            if value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            attributes[key] = value

    return {
        "gene_id": attributes.get("gene_id", "."),
        "gene_name": attributes.get("gene_name", "."),
        "gene_type": attributes.get("gene_type", "."),
    }


def parse_intersect_file(input_path: str) -> list[dict]:
    """
    Parse bedtools intersect output file.

    Expected format: 17 columns from `bedtools intersect -wa -wb -wo`
    - Columns 1-6: eccDNA BED6 (chr, start, end, name, score, strand)
    - Columns 7-16: Annotation BED (from gtf2bed)
    - Column 17: Overlap in base pairs

    Args:
        input_path: Path to the intersect output file.

    Returns:
        List of parsed row dictionaries.
    """
    parsed_data = []
    line_count = 0
    overlap_count = 0
    skipped_count = 0

    tab_splitter = re.compile(r"\t")

    with open(input_path) as f:
        for line in f:
            line_count += 1
            fields = tab_splitter.split(line.strip())

            if len(fields) != 17:
                if skipped_count < 20:
                    print(
                        f"Warning: Line {line_count} has {len(fields)} fields, "
                        f"expected 17. Skipping.",
                        file=sys.stderr,
                    )
                skipped_count += 1
                continue

            (ecc_chr, ecc_start_str, ecc_end_str,
             ecc_id, ecc_score, ecc_strand) = fields[0:6]

            (gene_chr, gene_start, gene_end, _, _, gene_strand,
             _, _, _, gene_attribute_string) = fields[6:16]

            overlap_bp_str = fields[16]

            # Check for valid overlap
            is_overlap = False
            if gene_chr != ".":
                try:
                    if float(overlap_bp_str) > 0:
                        is_overlap = True
                except ValueError:
                    pass

            if not is_overlap:
                continue

            overlap_count += 1

            gene_info = parse_attributes(gene_attribute_string)

            try:
                ecc_length = int(ecc_end_str) - int(ecc_start_str)
            except ValueError:
                ecc_length = "."

            row_data = {
                "ecc_chr": ecc_chr,
                "ecc_start": ecc_start_str,
                "ecc_end": ecc_end_str,
                "ecc_id": ecc_id,
                "ecc_length": ecc_length,
                "ecc_score": ecc_score,
                "ecc_strand": ecc_strand,
                "coverage": ".",
                "split_reads": ".",
                "discordant_mates": ".",
                "gene_id": gene_info["gene_id"],
                "gene_name": gene_info["gene_name"],
                "gene_type": gene_info["gene_type"],
                "gene_chr": gene_chr,
                "gene_start": gene_start,
                "gene_end": gene_end,
                "gene_strand": gene_strand,
                "overlap_bp": overlap_bp_str,
            }
            parsed_data.append(row_data)

    print(f"Processed {line_count} lines.")
    print(f"Found {overlap_count} overlapping entries.")
    if skipped_count > 0:
        print(f"Skipped {skipped_count} lines due to unexpected field count.")

    return parsed_data


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Parse bedtools intersect output for eccDNA-gene overlaps.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to bedtools intersect -wa -wb -wo output file",
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Path to output annotated TSV file",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Parsing: {args.input}")
    parsed_data = parse_intersect_file(args.input)

    if parsed_data:
        final_df = pd.DataFrame(parsed_data)
        final_df["ecc_start"] = pd.to_numeric(final_df["ecc_start"], errors="coerce")
        final_df["ecc_end"] = pd.to_numeric(final_df["ecc_end"], errors="coerce")
        final_df["ecc_score"] = pd.to_numeric(final_df["ecc_score"], errors="ignore")
        final_df["overlap_bp"] = pd.to_numeric(final_df["overlap_bp"], errors="coerce")

        final_df = final_df.reindex(columns=FINAL_COLUMNS, fill_value=".")
        final_df_sorted = final_df.sort_values(by=["ecc_chr", "ecc_start", "ecc_end"])
        final_df_sorted.to_csv(args.output, sep="\t", index=False, header=True, na_rep=".")

        print(f"\n✅ Created {args.output} ({len(final_df_sorted)} rows)")
    else:
        print("\nNo overlapping data found. Output file not created.")
        sys.exit(0)


if __name__ == "__main__":
    main()