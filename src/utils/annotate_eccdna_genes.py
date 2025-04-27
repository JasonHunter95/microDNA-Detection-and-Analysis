import pandas as pd
import re

def extract_gene_info(attribute_field):
    """
    Parse the GTF-style attribute column to extract gene_id, gene_name, and gene_type
    """
    gene_id = re.search(r'gene_id "([^"]+)"', attribute_field)
    gene_name = re.search(r'gene_name "([^"]+)"', attribute_field)
    gene_type = re.search(r'gene_type "([^"]+)"', attribute_field)

    return (
        gene_id.group(1) if gene_id else None,
        gene_name.group(1) if gene_name else None,
        gene_type.group(1) if gene_type else None
    )

def main():
    # load the intersected file
    df = pd.read_csv("data/SRR413969.eccdna.with_genes.bed", sep='\t', header=None)

    # parse relevant columns
    result = []
    for _, row in df.iterrows():
        ecc_id = row[3]
        chrom = row[0]
        start = row[1]
        end = row[2]
        attributes = row[15]

        gene_id, gene_name, gene_type = extract_gene_info(attributes)

        result.append([ecc_id, chrom, start, end, gene_id, gene_name, gene_type])

    # convert to DataFrame
    annotated = pd.DataFrame(result, columns=[
        "ecc_id", "chrom", "start", "end", "gene_id", "gene_name", "gene_type"
    ])

    # drop duplicates if any
    annotated = annotated.drop_duplicates()

    # save to TSV
    annotated.to_csv("data/SRR413969.eccdna.annotated.tsv", sep='\t', index=False)
    print("✅ Saved annotated file to data/SRR413969.eccdna.annotated.tsv")

if __name__ == "__main__":
    main()
