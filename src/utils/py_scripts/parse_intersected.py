import pandas as pd
import re
import sys
import os

# This needs to be parameterized in the future for every chromosome

intersect_file = '/Users/JasonHunter/Desktop/microDNA-Detection-and-Analysis/data/beds/chr1/SRR413984_chr1.eccdna.intersect_genes.bed'

output_dir = '/Users/JasonHunter/Desktop/microDNA-Detection-and-Analysis/data/results/chr1'
output_filename = 'SRR413984_chr1.intersected_annotations.tsv'
output_tsv = os.path.join(output_dir, output_filename)

os.makedirs(output_dir, exist_ok=True)

final_columns = [
    'ecc_chr', 'ecc_start', 'ecc_end', 'ecc_id', 'ecc_length',
    'ecc_score', 'ecc_strand', 'coverage', 'split_reads', 'discordant_mates', # Keep placeholders for now
    'gene_id', 'gene_name', 'gene_type',
    'gene_chr', 'gene_start', 'gene_end', 'gene_strand',
    'overlap_bp'
]

def parse_attributes(attr_string):
    attributes = {}
    if not isinstance(attr_string, str): return {'gene_id': '.', 'gene_name': '.', 'gene_type': '.'}
    for item in attr_string.strip().split(';'):
        item = item.strip();
        if not item: continue
        parts = item.split(' ', 1);
        if len(parts) == 2:
            key = parts[0].strip(); value = parts[1].strip()
            if value.startswith('"') and value.endswith('"'): value = value[1:-1]
            attributes[key] = value
    return { 'gene_id': attributes.get('gene_id', '.'), 'gene_name': attributes.get('gene_name', '.'), 'gene_type': attributes.get('gene_type', '.') }

parsed_data = []
line_count = 0
overlap_count = 0
skipped_count = 0

tab_splitter = re.compile(r'\t')

print(f"Attempting to read file: {intersect_file}")
if not os.path.exists(intersect_file):
    print(f"!!! ERROR: File not found at path: {intersect_file}")
    sys.exit(1)

try:
    with open(intersect_file, 'r') as f:
        for line in f:
            line_count += 1
            fields = tab_splitter.split(line.strip())

            if len(fields) != 17:
                 if skipped_count < 20:
                      print(f"Warning: Line {line_count} has {len(fields)} fields (split by regex), expected 17. Skipping: {line.strip()[:100]}...") # Print start of line
                 skipped_count += 1
                 continue

            ecc_chr, ecc_start_str, ecc_end_str, ecc_id, ecc_score, ecc_strand = fields[0:6]
            gene_chr, gene_start, gene_end, gene_name_field, gene_score, gene_strand, \
            gene_source, gene_feature, gene_frame, gene_attribute_string = fields[6:16]
            overlap_bp_str = fields[16]

            is_overlap = False
            if gene_chr != '.':
                try:
                    if float(overlap_bp_str) > 0:
                        is_overlap = True
                except ValueError:
                    pass # overlap_bp is not numeric

            if not is_overlap:
                continue

            overlap_count += 1

            gene_info = parse_attributes(gene_attribute_string)
            try:
                ecc_length = int(ecc_end_str) - int(ecc_start_str)
            except ValueError:
                ecc_length = '.'

            row_data = {
                'ecc_chr': ecc_chr, 'ecc_start': ecc_start_str, 'ecc_end': ecc_end_str, 'ecc_id': ecc_id, 'ecc_length': ecc_length,
                'ecc_score': ecc_score, 'ecc_strand': ecc_strand,
                'coverage': '.', 'split_reads': '.', 'discordant_mates': '.', # Add placeholders
                'gene_id': gene_info['gene_id'], 'gene_name': gene_info['gene_name'], 'gene_type': gene_info['gene_type'],
                'gene_chr': gene_chr, 'gene_start': gene_start, 'gene_end': gene_end, 'gene_strand': gene_strand,
                'overlap_bp': overlap_bp_str
            }
            parsed_data.append(row_data)

except Exception as e:
    print(f"An error occurred processing {intersect_file} around line {line_count}: {e}")
    raise e

print(f"\nProcessed {line_count} lines.")
print(f"Found {overlap_count} overlapping entries.")
if skipped_count > 0:
     print(f"Skipped {skipped_count} lines due to unexpected field count (expected 17).")


if parsed_data:
    final_df = pd.DataFrame(parsed_data)
    final_df['ecc_start'] = pd.to_numeric(final_df['ecc_start'], errors='coerce')
    final_df['ecc_end'] = pd.to_numeric(final_df['ecc_end'], errors='coerce')
    final_df['ecc_score'] = pd.to_numeric(final_df['ecc_score'], errors='ignore')
    final_df['overlap_bp'] = pd.to_numeric(final_df['overlap_bp'], errors='coerce')
    # reindex
    final_df = final_df.reindex(columns=final_columns, fill_value='.')
    final_df_sorted = final_df.sort_values(by=['ecc_chr', 'ecc_start', 'ecc_end'])
    final_df_sorted.to_csv(output_tsv, sep='\t', index=False, header=True, na_rep='.')
    print(f"\nSuccessfully created {output_tsv}")
else:
    print("\nNo overlapping data found after manual parsing. Output file not created.")