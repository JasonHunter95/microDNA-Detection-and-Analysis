# Filename: annotate_eccdna_closest.py
import pandas as pd
import re
import argparse # Import argparse
import os
import sys

def extract_gene_info(attribute_field):
    """
    Parse the GTF-style attribute column (typically from column 9+ in GTF,
    check which column it ends up in after gtf2bed conversion used as input B for closest)
    to extract multiple useful features. Robustly handles missing fields.
    """
    # Ensure input is treated as a string, even if it's None or NaN from pandas
    attribute_field = str(attribute_field) if attribute_field is not None else ""

    gene_id = re.search(r'gene_id "([^"]+)"', attribute_field)
    gene_name = re.search(r'gene_name "([^"]+)"', attribute_field)
    gene_type = re.search(r'gene_type "([^"]+)"', attribute_field)
    transcript_id = re.search(r'transcript_id "([^"]+)"', attribute_field)
    transcript_name = re.search(r'transcript_name "([^"]+)"', attribute_field)
    level = re.search(r'level (\d+)', attribute_field)
    gene_status = re.search(r'gene_status "([^"]+)"', attribute_field)
    # Add more extractions if needed

    return {
        "gene_id": gene_id.group(1) if gene_id else None,
        "gene_name": gene_name.group(1) if gene_name else None,
        "gene_type": gene_type.group(1) if gene_type else None,
        "transcript_id": transcript_id.group(1) if transcript_id else None,
        "transcript_name": transcript_name.group(1) if transcript_name else None,
        "level": int(level.group(1)) if level else None,
        "gene_status": gene_status.group(1) if gene_status else None
    }

def main():
    # --- Add Argument Parser ---
    parser = argparse.ArgumentParser(
        description="Annotate eccDNA based on bedtools closest -d results. Parses GTF attributes from the annotation file portion."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input file from 'bedtools closest -d'"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output annotated TSV file"
    )
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output

    print(f"Reading closest results from: {input_path}")
    print(f"Outputting annotated TSV to: {output_path}")

    try:
        # Define column names/indices based on bedtools closest -d output
        # Output: <cols from A> <cols from B> <distance>
        # Assumes input A (eccDNA) was BED6 (6 cols)
        # Assumes input B (annotation) was derived from GTF (e.g., via gtf2bed filtered for genes)
        # The number of cols from B depends on that conversion. Let's assume the original script's
        # indices were correct, implying file B had at least 10 relevant columns mapped to indices 6-15.
        # Col 16 (index 15) holds attributes, Col 17 (index 16) holds distance.

        # Load the file without a header
        df = pd.read_csv(input_path, sep='\t', header=None)

        # Check if the assumed attribute and distance columns exist
        attribute_col_index = 15 # Based on user's original code (row[15])
        distance_col_index = 16  # Based on user's original code (row[16]) and -d behavior

        if distance_col_index >= df.shape[1]:
             print(f"Error: Expected distance column at index {distance_col_index}, but file only has {df.shape[1]} columns.", file=sys.stderr)
             print("Did you run 'bedtools closest' with the '-d' option?", file=sys.stderr)
             sys.exit(1)
        if attribute_col_index >= df.shape[1]:
             print(f"Warning: Expected attribute column at index {attribute_col_index}, using None. File has {df.shape[1]} columns.", file=sys.stderr)
             # Allow processing to continue, but attributes will be None

        print(f"Parsing attributes from column index {attribute_col_index} (column number {attribute_col_index + 1}).")
        print(f"Parsing distance from column index {distance_col_index} (column number {distance_col_index + 1}).")


        # Parse relevant columns using indices from user's code
        result = []
        for _, row in df.iterrows():
            # Columns from eccDNA file (file A) - assuming BED6 input
            ecc_id = row[3]        # col 4: name
            ecc_chrom = row[0]     # col 1: chrom
            ecc_start = row[1]     # col 2: start
            ecc_end = row[2]       # col 3: end
            # Optional: ecc_score = row[4], ecc_strand = row[5]

            # Columns from Annotation file (file B) - indices based on user code
            # These might shift depending on exact BED used as input B for closest
            gene_chrom = row[6]    # col 7: anno_chrom
            gene_start = row[7]    # col 8: anno_start
            gene_end = row[8]      # col 9: anno_end
            gene_id_versioned = row[9] # col 10: anno_name (often complex)
            gene_strand = row[11]  # col 12: anno_strand (check if B file has strand)
            gene_source = row[12]  # col 13: anno_source (check if B file has source)
            gene_feature_type = row[13] # col 14: anno_feature (check if B file has feature)

            # Attribute column (check if exists)
            attributes = str(row[attribute_col_index]) if attribute_col_index < df.shape[1] else None

            # Distance column (already checked it exists)
            distance = row[distance_col_index]

            # Extract info from the attributes string
            extracted = extract_gene_info(attributes)

            # Append selected data
            result.append([
                ecc_id, ecc_chrom, ecc_start, ecc_end,
                # Optional: ecc_score, ecc_strand,
                gene_chrom, gene_start, gene_end,
                gene_id_versioned, # Keep the original name field from annotation BED
                extracted["gene_id"], extracted["gene_name"], extracted["gene_type"],
                extracted["transcript_id"], extracted["transcript_name"],
                extracted["level"], extracted["gene_status"],
                gene_strand, gene_source, gene_feature_type,
                distance # Add distance
            ])

        # Convert list of lists to DataFrame
        annotated_df = pd.DataFrame(result, columns=[
            "ecc_id", "ecc_chrom", "ecc_start", "ecc_end",
            # Optional: "ecc_score", "ecc_strand",
            "gene_chrom", "gene_start", "gene_end",
            "gene_identifier_anno", # Renamed original col 10
            "gene_id", "gene_name", "gene_type",
            "transcript_id", "transcript_name",
            "level", "gene_status",
            "gene_strand", "gene_source", "gene_feature_type",
            "distance" # Add distance column name
        ])

        # Drop duplicate rows (if an eccDNA is equidistant from the same annotation feature reported multiple times)
        annotated_df = annotated_df.drop_duplicates()

        # save to TSV using the parameterized output path
        annotated_df.to_csv(output_path, sep='\t', index=False)

        print(f"✅ Saved annotated file ({len(annotated_df)} rows) to: {os.path.abspath(output_path)}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Warning: Input file {input_path} is empty. No output generated.", file=sys.stderr)
        # Create an empty output file with header
        header_cols = [
            "ecc_id", "ecc_chrom", "ecc_start", "ecc_end",
            "gene_chrom", "gene_start", "gene_end", "gene_identifier_anno",
            "gene_id", "gene_name", "gene_type", "transcript_id", "transcript_name",
            "level", "gene_status", "gene_strand", "gene_source", "gene_feature_type", "distance"
        ]
        try:
            with open(output_path, 'w') as outfile:
               outfile.write('\t'.join(header_cols) + '\n')
            print(f"Created empty output file with header: {os.path.abspath(output_path)}")
        except Exception as e:
            print(f"Error creating empty output file {output_path}: {e}", file=sys.stderr)
        sys.exit(0) # Successful exit, just no data
    except Exception as e:
        print(f"An error occurred during annotation processing: {e}", file=sys.stderr)
        # Attempt to clean up potentially incomplete output file
        if os.path.exists(output_path):
             try:
                 os.remove(output_path)
                 print(f"Removed potentially incomplete output file: {output_path}", file=sys.stderr)
             except OSError as oe:
                 print(f"Error trying to remove incomplete output file {output_path}: {oe}", file=sys.stderr)
        sys.exit(1) # Exit with error

if __name__ == "__main__":
    main()