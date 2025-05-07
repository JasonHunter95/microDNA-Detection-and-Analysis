import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_distribution(data, column, title, xlabel, ylabel):
    """
    Plot the distribution of a given column in the DataFrame.
    
    Parameters:
    - data: DataFrame containing the data
    - column: Column name to plot
    - title: Title of the plot
    - xlabel: Label for the x-axis
    - ylabel: Label for the y-axis
    """
    fig = plt.figure(figsize=(12, 6))
    
    data_to_plot = data[data['ecc_score'] > 0]
    sns.histplot(data=data_to_plot, x='ecc_score', bins='auto', kde=True, log_scale=True) # Use log_scale=True
    plt.title('Distribution of eccDNA Confidence Scores (Log Scale)')
    plt.xlabel('Score [Log Scale]')
    plt.ylabel('Frequency')
    
    plt.axvline(data[column].mean(), color='r', linestyle='--', label='Mean')
    plt.axvline(data[column].median(), color='g', linestyle='--', label='Median')
    plt.axvline(data[column].quantile(0.75), color='b', linestyle='--', label='75th Percentile')
    plt.axvline(data[column].quantile(0.25), color='y', linestyle='--', label='25th Percentile')
    plt.legend()
    return fig
    
def main():
    data_path = "data/results/chr1/SRR413984_cleaner_chr1.unique_eccdna_table.tsv"
    ecc_data = pd.read_csv(data_path, sep='\t', low_memory=False)
    
    if 'ecc_score' not in ecc_data.columns:
        print("Required column 'ecc_score' missing from the data.")
        return
    
    plot_distribution(ecc_data, 'ecc_score', 'Distribution of eccDNA Scores', 'Score', 'Frequency')
    
    output_path = 'figures/'
    plt.savefig(output_path + 'eccDNA_score_distribution.png')
    plt.show()
    plt.close()
    
if __name__ == "__main__":
    main()