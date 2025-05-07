import pandas as pd

df = pd.read_csv("data/results/chr1/SRR413984_chr1.intersected_annotations.tsv", sep="\t", header=0)
# drop the empty columns
df = df.drop(columns=["ecc_strand","coverage","split_reads","discordant_mates"])
# write back out
df.to_csv("data/results/chr1/SRR413984_cleaner_chr1.intersected_annotations.tsv", sep="\t", index=False)
