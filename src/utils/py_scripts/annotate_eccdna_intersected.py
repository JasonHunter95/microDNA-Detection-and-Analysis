# Filename: annotate_eccdna_intersected.py (Assumed name)
import pandas as pd
import re
import argparse # Import argparse
import os
import sys

def extract_gene_info(attribute_field):
    """
    Parse the GTF-style attribute column (typically from column 9+ in GTF,
    so check which column it ends up in after bedtools intersect -wb)
    to extract multiple useful features. Robustly handles missing fields.
    """
    gene_id = re.search(r'gene_id "([^"]+)"', attribute_field)
    gene_name = re.search(r'gene_name "([^"]+)"', attribute_field)
    gene_type = re.search(r'gene_type "([^"]+)"', attribute_field)
    transcript_id = re.search(r'transcript_id "([^"]+)"', attribute_field)
    transcript_name = re.search(r'transcript_name "([^"]+)"', attribute_field)
    level = re.search(r'level (\d+)', attribute_field)
    gene_status = re.search(r'gene_status "([^"]+)"', attribute_field)
    # Add more extractions if needed (e.g., exon_id, exon_number)

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
        description="Annotate eccDNA based on bedtools intersect results (-wa -wb). Parses GTF attributes from the annotation file portion."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input file from 'bedtools intersect -wa -wb'"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output annotated TSV file"
    )
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output

    print(f"Reading intersect results from: {input_path}")
    print(f"Outputting annotated TSV to: {output_path}")

    try:
        # Define column names based on bedtools intersect -wa -wb
        # Assumes input A (eccDNA) was BED6 and input B (annotation) was BED-like from gtf2bed
        # Adjust num_a_cols if your input A has a different number of columns
        num_a_cols = 6 # e.g., chr, start, end, name, score, strand
        # The columns from B will depend on gtf2bed output. GTF has 9 columns + attributes.
        # gtf2bed might simplify this. Let's assume the user's original parsing was correct
        # and the important attributes are in the 16th column (index 15) overall.

        # Load the file without a header
        df = pd.read_csv(input_path, sep='\t', header=None)

        # Check if the assumed attribute column index exists
        attribute_col_index = 15 # Based on user's original code (row[15])
        if attribute_col_index >= df.shape[1]:
             print(f"Error: Expected attribute column at index {attribute_col_index}, but file only has {df.shape[1]} columns.", file=sys.stderr)
             print("Check the output of your 'bedtools intersect -wa -wb' command.", file=sys.stderr)
             sys.exit(1)

        print(f"Parsing attributes from column index {attribute_col_index} (column number {attribute_col_index + 1}).")

        # parse relevant columns using indices from user's code
        result = []
        for _, row in df.iterrows():
            # Columns from eccDNA file (file A) - assuming BED6 input
            ecc_id = row[3]        # col 4: name
            ecc_chrom = row[0]     # col 1: chrom
            ecc_start = row[1]     # col 2: start
            ecc_end = row[2]       # col 3: end
            # Optional: ecc_score = row[4], ecc_strand = row[5]

            # Columns from Annotation file (file B) - indices based on user code
            # These might shift depending on exact gtf2bed output used with intersect
            gene_chrom = row[6]    # col 7: anno_chrom
            gene_start = row[7]    # col 8: anno_start
            gene_end = row[8]      # col 9: anno_end
            gene_id_versioned = row[9] # col 10: anno_name (often complex)
            gene_strand = row[11]  # col 12: anno_strand
            gene_source = row[12]  # col 13: anno_source (from GTF col 2)
            gene_feature_type = row[13] # col 14: anno_feature (from GTF col 3)
            attributes = str(row[attribute_col_index]) # col 16: anno_attributes (GTF col 9)

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
                gene_strand, gene_source, gene_feature_type
            ])

        # convert list of lists to DataFrame
        annotated_df = pd.DataFrame(result, columns=[
            "ecc_id", "ecc_chrom", "ecc_start", "ecc_end",
            # Optional: "ecc_score", "ecc_strand",
            "gene_chrom", "gene_start", "gene_end",
            "gene_identifier_anno", # Renamed original col 10
            "gene_id", "gene_name", "gene_type",
            "transcript_id", "transcript_name",
            "level", "gene_status",
            "gene_strand", "gene_source", "gene_feature_type"
        ])

        # drop duplicate rows (if an eccDNA overlaps the same annotation feature multiple times in the input)
        annotated_df = annotated_df.drop_duplicates()

        # save to TSV using the parameterized output path
        annotated_df.to_csv(output_path, sep='\t', index=False)

        print(f"✅ Saved annotated file ({len(annotated_df)} rows) to: {os.path.abspath(output_path)}")

    except FileNotFoundError:
        print(f"Error: Input file not found at {input_path}", file=sys.stderr)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Warning: Input file {input_path} is empty. No output generated.", file=sys.stderr)
        # Create an empty output file with header perhaps?
        # Or just exit gracefully. For now, just print warning.
        # To create empty file with header:
        # header_cols = ["ecc_id", "ecc_chrom", ...] # list all cols
        # with open(output_path, 'w') as outfile:
        #    outfile.write('\t'.join(header_cols) + '\n')
        # print(f"Created empty output file with header: {os.path.abspath(output_path)}")
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