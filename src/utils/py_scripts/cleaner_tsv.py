import pandas as pd
dtype_spec = {
    'ecc_start': 'int64',
    'ecc_end': 'int64',
    'ecc_length': 'int64',
    'ecc_score': 'float64',
    'gene_start': 'int64',
    'gene_end': 'int64',
    'overlap_bp': 'int64',
    'ecc_chr': 'category',
    'gene_chr': 'category',
    'gene_strand': 'category',
    'gene_type': 'category',
}

chunk_size = 100000
chunks = []
try:
    for chunk in pd.read_csv('data/results/chr1/SRR413984_cleaner_chr1.intersected_annotations.tsv', sep='\t', dtype=dtype_spec, chunksize=chunk_size, low_memory=False):
        chunks.append(chunk)
    df = pd.concat(chunks, ignore_index=True)
    print("File loaded successfully.")
except FileNotFoundError:
    print("Error: File not found. Make sure the path is correct.")
    exit()
except Exception as e:
    print(f"An error occurred during loading: {e}")
    try:
        print("Attempting load without specified dtypes...")
        chunks = []
        for chunk in pd.read_csv('data/results/chr1/SRR413984_cleaner_chr1.intersected_annotations.tsv', sep='\t', chunksize=chunk_size, low_memory=False):
             chunks.append(chunk)
        df = pd.concat(chunks, ignore_index=True)
        print("File loaded successfully (without initial dtypes).")
    except Exception as e2:
        print(f"Loading failed again: {e2}")
        exit()

print(f"Initial DataFrame shape: {df.shape}")
print(df.head())
print(df.info())
numeric_cols = ['ecc_start', 'ecc_end', 'ecc_length', 'ecc_score', 'gene_start', 'gene_end', 'overlap_bp']
for col in numeric_cols:
    if df[col].dtype == 'object':
        df[col] = pd.to_numeric(df[col], errors='coerce')
        
print("\nMissing values per column:")
print(df.isnull().sum())

critical_cols = ['ecc_chr', 'ecc_start', 'ecc_end', 'ecc_id', 'gene_id', 'overlap_bp']
df.dropna(subset=critical_cols, inplace=True)

initial_rows = len(df)
df.drop_duplicates(inplace=True)
print(f"\nRemoved {initial_rows - len(df)} duplicate rows.")

print("\nValidating coordinates...")
invalid_ecc_coords = df[df['ecc_start'] > df['ecc_end']].shape[0]
invalid_gene_coords = df[df['gene_start'] > df['gene_end']].shape[0]

expected_length = df['ecc_end'] - df['ecc_start'] + 1
length_diff = abs(df['ecc_length'] - expected_length)
inconsistent_length = df[length_diff > 1].shape[0]

gene_length = df['gene_end'] - df['gene_start'] + 1
invalid_overlap_ecc = df[df['overlap_bp'] > df['ecc_length']].shape[0]

valid_gene_rows = df.dropna(subset=['gene_start', 'gene_end'])
invalid_overlap_gene = valid_gene_rows[valid_gene_rows['overlap_bp'] > (valid_gene_rows['gene_end'] - valid_gene_rows['gene_start'] + 1)].shape[0]


print(f"Rows with invalid ecc coords (start > end): {invalid_ecc_coords}")
print(f"Rows with invalid gene coords (start > end): {invalid_gene_coords}")
print(f"Rows with inconsistent ecc_length (diff > 1bp): {inconsistent_length}")
print(f"Rows where overlap_bp > ecc_length: {invalid_overlap_ecc}")
print(f"Rows where overlap_bp > gene_length: {invalid_overlap_gene}")


print(f"\necc_score distribution:\n{df['ecc_score'].describe()}")

score_threshold = 10
df_filtered = df[df['ecc_score'] >= score_threshold].copy()
print(f"\nRows remaining after filtering ecc_score >= {score_threshold}: {len(df_filtered)}")

df_filtered['gene_length'] = df_filtered['gene_end'] - df_filtered['gene_start'] + 1
df_filtered['perc_overlap_ecc'] = df_filtered.apply(
    lambda row: (row['overlap_bp'] / row['ecc_length'] * 100) if row['ecc_length'] > 0 else 0, axis=1
)
df_filtered['perc_overlap_gene'] = df_filtered.apply(
    lambda row: (row['overlap_bp'] / row['gene_length'] * 100) if pd.notna(row['gene_length']) and row['gene_length'] > 0 else 0, axis=1
)

min_bp_overlap = 50
min_perc_overlap_ecc = 10

df_meaningful_overlap = df_filtered[
    (df_filtered['overlap_bp'] >= min_bp_overlap) &
    (df_filtered['perc_overlap_ecc'] >= min_perc_overlap_ecc)
].copy()
print(f"\nRows remaining after filtering for meaningful overlap: {len(df_meaningful_overlap)}")

output_path = 'data/results/chr1/SRR413984_cleaner_chr1.meaningful_data.tsv'
df_meaningful_overlap.to_csv(output_path, sep='\t', index=False)
print(f"\nCleaned data saved to {output_path}")