#!/usr/bin/env python3
"""
clean_bed.py - Format BED-like files into BED6 format with unique eccDNA names.

Processes Circle-Map output or other BED-like files to produce standard BED6
format with unique identifiers for each eccDNA candidate.
"""

import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Format a BED-like file into BED6 format with unique eccDNA names. Processes line-by-line for low memory usage."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input BED-like file (e.g., Circle-Map output or cleaned BED)" # Modified help slightly
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output BED6 formatted file"
    )
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output
    ecc_counter = 1 # Counter for generating unique names

    print(f"Processing file: {input_path}")
    print(f"Outputting BED6 to: {output_path}")

    try:
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line_num, line in enumerate(infile, 1):
                # Skip header or comment lines
                if line.startswith('#'):
                    continue

                # Skip empty lines
                stripped_line = line.strip()
                if not stripped_line:
                    continue

                fields = stripped_line.split('\t')

                # --- Input Flexibility: Adapt based on expected input ---
                # Expecting at least chrom, start, end from cleaned Circle-Map BED
                # Or potentially more if input is from intersect/closest
                if len(fields) < 3:
                    print(f"Warning: Skipping line {line_num}. Expected at least 3 columns (chr, start, end), found {len(fields)}: {stripped_line}", file=sys.stderr)
                    continue

                # --- Extract core BED fields ---
                chrom = fields[0]
                start = fields[1]
                end = fields[2]

                # --- Assign Score ---
                # Use column 4 if it exists and looks numeric-ish, otherwise default to '0'
                score = '0' # Default score
                if len(fields) > 3:
                     # Basic check if it could be a score (integer or float representation)
                     if fields[3].replace('.', '', 1).isdigit() or (fields[3].startswith('-') and fields[3][1:].replace('.', '', 1).isdigit()):
                          score = fields[3]
                     # Optional: If col 4 might be the name from an intermediate step,
                     # you might check `len(fields) > 4` for the score instead.

                # --- Assign Strand ---
                # Use column 6 if it exists and is '+' or '-', otherwise default to '.'
                strand = '.' # Default strand
                if len(fields) > 5 and fields[5] in ['+', '-']:
                    strand = fields[5]

                # Generate the unique eccDNA name
                name = f"ecc_{str(ecc_counter).zfill(5)}"

                # Format the output BED6 line: chr, start, end, name, score, strand
                output_line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"

                # Write to the output file
                outfile.write(output_line)

                ecc_counter += 1

        processed_count = ecc_counter - 1
        print(f"✅ Successfully formatted {processed_count} records.")
        print(f"✅ Saved BED6 file to: {os.path.abspath(output_path)}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}", file=sys.stderr)
        sys.exit(1) # Exit with a non-zero status to indicate error
    except Exception as e:
        print(f"An error occurred during processing: {e}", file=sys.stderr)
        # Attempt to clean up potentially incomplete output file upon error
        if os.path.exists(output_path):
             try:
                 os.remove(output_path)
                 print(f"Removed potentially incomplete output file: {output_path}", file=sys.stderr)
             except OSError as oe:
                 print(f"Error trying to remove incomplete output file {output_path}: {oe}", file=sys.stderr)
        sys.exit(1) # Exit with a non-zero status

if __name__ == "__main__":
    main()