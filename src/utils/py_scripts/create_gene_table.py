import pandas as pd

gene_table = pd.read_csv("data/results/chr1/SRR413984_cleaner_chr1.meaningful_data.tsv", sep='\t', low_memory=False)

print("Available columns:")
print(gene_table.columns.tolist())

print("\nFirst few rows:")
print(gene_table.head())

gene_centric = gene_table.dropna(subset=['gene_id']).groupby('gene_id').agg(
    gene_name=('gene_name', 'first'),
    gene_type=('gene_type', 'first'),
    gene_chr=('gene_chr', 'first'),
    gene_start=('gene_start', 'first'),
    gene_end=('gene_end', 'first'),
    gene_strand=('gene_strand', 'first'),
    num_overlapping_eccdnas=('ecc_id', 'nunique'),
    overlapping_ecc_ids=('ecc_id', lambda x: list(x.unique())),
    avg_ecc_score=('ecc_score', 'mean'),
    median_ecc_length=('ecc_length', 'median'),
    total_overlap_bp=('overlap_bp', 'sum')
).reset_index()

print("\nGene-centric table sample:")
print(gene_centric.sort_values('num_overlapping_eccdnas', ascending=False).head())
print(f"Shape of gene-centric table: {gene_centric.shape}")
output_path = 'data/results/chr1/SRR413984_cleaner_chr1.gene_centric_table.tsv'
gene_centric.to_csv(output_path, sep='\t', index=False)
print(f"\nGene-centric table saved to {output_path}")