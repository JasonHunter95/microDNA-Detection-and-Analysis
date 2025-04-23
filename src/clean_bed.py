import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Format raw Circle-Map output (BED-like) into BED6 with eccDNA names.")
    parser.add_argument("-i", "--input", required=True, help="Path to raw Circle-Map output file (e.g., peaks.bed)")
    parser.add_argument("-o", "--output", required=True, help="Output BED6 file")
    args = parser.parse_args()

    # read file as tab-separated with no header
    try:
        df = pd.read_csv(args.input, sep='\t', comment='#', header=None)
    except Exception as e:
        raise ValueError(f"Failed to read input file: {e}")

    # check minimum number of columns
    if df.shape[1] < 3:
        raise ValueError("Input must have at least 3 columns: chrom, start, end")

    # extract first 4 columns and name/strand placeholders
    df = df.iloc[:, :4].copy()
    df.columns = ['chrom', 'start', 'end', 'score']
    df['name'] = [f"ecc_{i:05d}" for i in range(1, len(df) + 1)]
    df['strand'] = '+'

    # reorder to BED6 and write output
    df = df[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    df.to_csv(args.output, sep='\t', header=False, index=False)
    print(f"✅ Saved cleaned BED6 file to: {os.path.abspath(args.output)}")

if __name__ == "__main__":
    main()
