import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("data/SRR413969.eccdna.annotated.tsv", sep='\t')
df['length'] = df['end'] - df['start']

plt.hist(df['length'], bins=50)
plt.xlabel("eccDNA length (bp)")
plt.ylabel("Count")
plt.title("eccDNA Size Distribution")
plt.show()
plt.savefig("figures/eccdna_size_distribution.png")