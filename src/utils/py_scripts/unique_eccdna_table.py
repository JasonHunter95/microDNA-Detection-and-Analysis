import pandas as pd

df_eccdna_unique = pd.read_csv("data/results/chr1/SRR413984_cleaner_chr1.intersected_annotations.tsv", sep="\t", header=0)

eccdna_unique = df_eccdna_unique.groupby('ecc_id').agg(
    ecc_chr=('ecc_chr', 'first'),
    ecc_start=('ecc_start', 'first'),
    ecc_end=('ecc_end', 'first'),
    ecc_length=('ecc_length', 'first'),
    ecc_score=('ecc_score', 'first'),
    num_overlapping_genes=('gene_id', 'nunique'), # Count unique genes overlapping this eccDNA
    overlapping_gene_ids=('gene_id', lambda x: list(x.dropna().unique())), # List unique gene IDs
    overlapping_gene_names=('gene_name', lambda x: list(x.dropna().unique())), # List unique gene names
    overlapping_gene_types=('gene_type', lambda x: list(x.dropna().unique())) # List unique gene types
).reset_index()

print("\nUnique eccDNA table sample:")
print(eccdna_unique.head())
print(f"Shape of unique eccDNA table: {eccdna_unique.shape}")

output_path = 'data/results/chr1/SRR413984_cleaner_chr1.unique_eccdna_table.tsv'
eccdna_unique.to_csv(output_path, sep='\t', index=False)
print(f"\nUnique eccDNA table saved to {output_path}")