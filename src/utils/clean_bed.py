import argparse
import os
import sys # For writing error messages to stderr

def main():
    parser = argparse.ArgumentParser(
        description="Format a BED-like file into BED6 format with unique eccDNA names. Processes line-by-line for low memory usage."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input BED-like file (e.g., Circle-Map output or intersect results)"
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
    print(f"Outputting to: {output_path}")

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

                # Ensure we have at least the first 4 columns needed (chrom, start, end, score/placeholder)
                # Adjust this number if your input format guarantees fewer/more columns initially
                if len(fields) < 4:
                    print(f"Warning: Skipping line {line_num}. Expected at least 4 columns, found {len(fields)}: {stripped_line}", file=sys.stderr)
                    continue

                # Extract required fields
                chrom = fields[0]
                start = fields[1]
                end = fields[2]
                # Use the 4th column as the score. If it doesn't exist or isn't numeric,
                # BED format often uses '0' or '1000'. Let's use the value directly if present.
                score = fields[3]
                # BED6 requires a strand column. Use '.' if not present or if input isn't stranded.
                strand = '.'

                # Generate the unique eccDNA name
                # Using zfill for consistent padding (e.g., ecc_00001)
                name = f"ecc_{str(ecc_counter).zfill(5)}"

                # Format the output BED6 line
                # Ensure fields are tab-separated
                output_line = f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"

                # Write to the output file
                outfile.write(output_line)

                ecc_counter += 1

        processed_count = ecc_counter - 1
        print(f"✅ Successfully processed {processed_count} records.")
        print(f"✅ Saved cleaned BED6 file to: {os.path.abspath(output_path)}")

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