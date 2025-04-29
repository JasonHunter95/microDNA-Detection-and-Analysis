import pandas as pd
import re

def extract_gene_info(attribute_field):
    """
    Parse the GTF-style attribute column to extract multiple useful features
    """
    gene_id = re.search(r'gene_id "([^"]+)"', attribute_field)
    gene_name = re.search(r'gene_name "([^"]+)"', attribute_field)
    gene_type = re.search(r'gene_type "([^"]+)"', attribute_field)
    transcript_id = re.search(r'transcript_id "([^"]+)"', attribute_field)
    transcript_name = re.search(r'transcript_name "([^"]+)"', attribute_field)
    level = re.search(r'level (\d+)', attribute_field)
    gene_status = re.search(r'gene_status "([^"]+)"', attribute_field)

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
    # load the bedtools closest file
    df = pd.read_csv("data/beds/chr1/SRR413984_chr1.eccdna.intersected_genes.bed", sep='\t', header=None)

    # parse relevant columns
    result = []
    for _, row in df.iterrows():
        ecc_id = row[3]
        ecc_chrom = row[0]
        ecc_start = row[1]
        ecc_end = row[2]
        gene_chrom = row[6]
        gene_start = row[7]
        gene_end = row[8]
        gene_id_versioned = row[9]
        gene_strand = row[11]
        gene_source = row[12]
        gene_feature_type = row[13]
        attributes = row[15]

        extracted = extract_gene_info(attributes)

        result.append([
            ecc_id, ecc_chrom, ecc_start, ecc_end,
            gene_chrom, gene_start, gene_end,
            gene_id_versioned,
            extracted["gene_id"], extracted["gene_name"], extracted["gene_type"],
            extracted["transcript_id"], extracted["transcript_name"],
            extracted["level"], extracted["gene_status"],
            gene_strand, gene_source, gene_feature_type
        ])

    # convert to DataFrame
    annotated = pd.DataFrame(result, columns=[
        "ecc_id", "ecc_chrom", "ecc_start", "ecc_end",
        "gene_chrom", "gene_start", "gene_end",
        "gene_id_versioned",
        "gene_id", "gene_name", "gene_type",
        "transcript_id", "transcript_name",
        "level", "gene_status",
        "gene_strand", "gene_source", "gene_feature_type"
    ])

    # drop duplicates if any
    annotated = annotated.drop_duplicates()

    # save to TSV
    annotated.to_csv("data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv", sep='\t', index=False)
    print("✅ Saved annotated file to data/results/chr1/SRR413984_chr1_intersected.eccdna.annotated.tsv")

if __name__ == "__main__":
    main()