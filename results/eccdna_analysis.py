#!/usr/bin/env python3
"""
Extended eccDNA Analysis Script

This script provides comprehensive analysis of eccDNA detection results including:
- Detailed size distribution statistics
- microDNA vs eccDNA characterization
- Gene type enrichment analysis with statistical testing
- Hotspot identification
- Distance-to-gene analysis
- Publication-quality figures

Usage:
    python eccdna_analysis.py <annotated_tsv> [--output-dir <dir>]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns

# Set publication-quality plot style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'figure.figsize': (10, 6),
    'figure.dpi': 150,
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'savefig.bbox': 'tight',
    'savefig.dpi': 300
})


def load_data(filepath: str) -> pd.DataFrame:
    """Load annotated eccDNA TSV file."""
    df = pd.read_csv(filepath, sep='\t')
    # Calculate size for each unique eccDNA
    df['size'] = df['ecc_end'] - df['ecc_start']
    return df


def get_unique_eccdna(df: pd.DataFrame) -> pd.DataFrame:
    """Get unique eccDNA entries (one row per ecc_id)."""
    return df.drop_duplicates('ecc_id')


def analyze_size_distribution(unique_ecc: pd.DataFrame, output_dir: Path) -> dict:
    """
    Comprehensive size distribution analysis.
    
    Returns statistics dictionary.
    """
    sizes = unique_ecc['size']
    
    # Calculate statistics
    stats_dict = {
        'count': len(sizes),
        'mean': sizes.mean(),
        'median': sizes.median(),
        'std': sizes.std(),
        'min': sizes.min(),
        'max': sizes.max(),
        'q25': sizes.quantile(0.25),
        'q75': sizes.quantile(0.75),
        'microDNA_count': (sizes < 400).sum(),
        'eccDNA_count': (sizes >= 400).sum(),
        'microDNA_pct': (sizes < 400).mean() * 100,
    }
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. Histogram with log scale
    ax1 = axes[0, 0]
    ax1.hist(sizes, bins=100, log=True, color='steelblue', edgecolor='white', alpha=0.8)
    ax1.axvline(400, color='red', linestyle='--', linewidth=2, label='microDNA threshold (400bp)')
    ax1.axvline(stats_dict['median'], color='orange', linestyle='-', linewidth=2, label=f'Median ({stats_dict["median"]:.0f}bp)')
    ax1.set_xlabel('eccDNA Size (bp)')
    ax1.set_ylabel('Count (log scale)')
    ax1.set_title('eccDNA Size Distribution')
    ax1.legend()
    
    # 2. Box plot comparing microDNA vs eccDNA
    ax2 = axes[0, 1]
    micro = sizes[sizes < 400]
    ecc = sizes[sizes >= 400]
    bp = ax2.boxplot([micro, ecc], labels=['microDNA\n(<400bp)', 'eccDNA\n(≥400bp)'], 
                     patch_artist=True)
    bp['boxes'][0].set_facecolor('lightcoral')
    bp['boxes'][1].set_facecolor('steelblue')
    ax2.set_ylabel('Size (bp)')
    ax2.set_title('Size Comparison: microDNA vs eccDNA')
    
    # Add counts as text
    ax2.text(1, micro.max() * 0.9, f'n={len(micro):,}', ha='center', fontsize=11)
    ax2.text(2, ecc.max() * 0.9, f'n={len(ecc):,}', ha='center', fontsize=11)
    
    # 3. Cumulative distribution
    ax3 = axes[1, 0]
    sorted_sizes = np.sort(sizes)
    cumulative = np.arange(1, len(sorted_sizes) + 1) / len(sorted_sizes)
    ax3.plot(sorted_sizes, cumulative, color='steelblue', linewidth=2)
    ax3.axvline(400, color='red', linestyle='--', alpha=0.7, label='400bp threshold')
    ax3.axhline(0.5, color='gray', linestyle=':', alpha=0.5)
    ax3.set_xlabel('eccDNA Size (bp)')
    ax3.set_ylabel('Cumulative Proportion')
    ax3.set_title('Cumulative Size Distribution')
    ax3.set_xscale('log')
    ax3.legend()
    
    # 4. Pie chart with improved aesthetics
    ax4 = axes[1, 1]
    labels = ['microDNA (<400bp)', 'eccDNA (≥400bp)']
    sizes_pie = [stats_dict['microDNA_count'], stats_dict['eccDNA_count']]
    colors = ['lightcoral', 'steelblue']
    explode = (0.02, 0)
    wedges, texts, autotexts = ax4.pie(sizes_pie, labels=labels, autopct='%1.1f%%',
                                        colors=colors, explode=explode, startangle=90,
                                        textprops={'fontsize': 11})
    ax4.set_title('Proportion of microDNA vs eccDNA')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'size_distribution_comprehensive.png')
    plt.close()
    
    print("\n" + "="*60)
    print("SIZE DISTRIBUTION STATISTICS")
    print("="*60)
    print(f"Total eccDNA:     {stats_dict['count']:,}")
    print(f"Mean size:        {stats_dict['mean']:,.1f} bp")
    print(f"Median size:      {stats_dict['median']:,.1f} bp")
    print(f"Std deviation:    {stats_dict['std']:,.1f} bp")
    print(f"Range:            {stats_dict['min']:,} - {stats_dict['max']:,} bp")
    print(f"IQR:              {stats_dict['q25']:,.1f} - {stats_dict['q75']:,.1f} bp")
    print(f"\nmicroDNA (<400bp): {stats_dict['microDNA_count']:,} ({stats_dict['microDNA_pct']:.1f}%)")
    print(f"eccDNA (≥400bp):   {stats_dict['eccDNA_count']:,} ({100-stats_dict['microDNA_pct']:.1f}%)")
    
    return stats_dict


def analyze_gene_type_enrichment(df: pd.DataFrame, output_dir: Path) -> pd.DataFrame:
    """
    Perform gene type enrichment analysis with Chi-squared test.
    
    Compares observed gene type frequencies in eccDNA vs expected background.
    """
    # Get gene type counts for unique eccDNA
    unique_ecc = get_unique_eccdna(df)
    gene_type_counts = unique_ecc['gene_type'].value_counts()
    
    # Calculate proportions
    total = gene_type_counts.sum()
    proportions = gene_type_counts / total
    
    # Create enrichment dataframe
    enrichment_df = pd.DataFrame({
        'gene_type': gene_type_counts.index,
        'count': gene_type_counts.values,
        'proportion': proportions.values
    })
    
    # For statistical testing, we'd need background genome proportions
    # Here we test against uniform distribution as a demonstration
    # In practice, you'd compare against the actual genomic proportions
    
    # Chi-squared test against uniform distribution (top 10 types)
    top_10 = gene_type_counts.head(10)
    expected = np.full(len(top_10), top_10.sum() / len(top_10))
    chi2, p_value = stats.chisquare(top_10.values, expected)
    
    print("\n" + "="*60)
    print("GENE TYPE ENRICHMENT ANALYSIS")
    print("="*60)
    print("\nTop 10 Gene Types Associated with eccDNA:")
    print("-" * 40)
    for gene_type, count in top_10.items():
        pct = count / total * 100
        print(f"  {gene_type:30s} {count:>6,} ({pct:>5.1f}%)")
    
    print(f"\nChi-squared test (vs uniform, top 10):")
    print(f"  χ² = {chi2:.2f}, p-value = {p_value:.2e}")
    if p_value < 0.001:
        print("  → Highly significant non-uniform distribution")
    
    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 1. Horizontal bar chart of top gene types
    ax1 = axes[0]
    top_15 = gene_type_counts.head(15)
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(top_15)))
    bars = ax1.barh(range(len(top_15)), top_15.values, color=colors)
    ax1.set_yticks(range(len(top_15)))
    ax1.set_yticklabels(top_15.index)
    ax1.set_xlabel('Count')
    ax1.set_title('Top 15 Gene Types Associated with eccDNA')
    ax1.invert_yaxis()
    
    # Add count labels
    for i, (count, bar) in enumerate(zip(top_15.values, bars)):
        ax1.text(count + max(top_15) * 0.01, i, f'{count:,}', va='center', fontsize=9)
    
    # 2. Stacked bar comparing microDNA vs eccDNA gene type preferences
    ax2 = axes[1]
    unique_ecc['is_microDNA'] = unique_ecc['size'] < 400
    
    # Get top 8 gene types for comparison
    top_8_types = gene_type_counts.head(8).index
    micro_counts = unique_ecc[unique_ecc['is_microDNA']]['gene_type'].value_counts()
    ecc_counts = unique_ecc[~unique_ecc['is_microDNA']]['gene_type'].value_counts()
    
    # Normalize to proportions
    micro_prop = micro_counts.reindex(top_8_types, fill_value=0) / micro_counts.sum() * 100
    ecc_prop = ecc_counts.reindex(top_8_types, fill_value=0) / ecc_counts.sum() * 100
    
    x = np.arange(len(top_8_types))
    width = 0.35
    
    ax2.bar(x - width/2, micro_prop.values, width, label='microDNA (<400bp)', color='lightcoral')
    ax2.bar(x + width/2, ecc_prop.values, width, label='eccDNA (≥400bp)', color='steelblue')
    
    ax2.set_xticks(x)
    ax2.set_xticklabels(top_8_types, rotation=45, ha='right')
    ax2.set_ylabel('Proportion (%)')
    ax2.set_title('Gene Type Distribution: microDNA vs eccDNA')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'gene_type_enrichment.png')
    plt.close()
    
    return enrichment_df


def identify_hotspots(df: pd.DataFrame, output_dir: Path, top_n: int = 20) -> pd.DataFrame:
    """
    Identify eccDNA hotspots - genes with the most overlapping eccDNA.
    """
    # Filter for genic eccDNA (distance == 0)
    genic = df[df['distance'] == 0]
    
    # Count eccDNA per gene
    gene_counts = genic.groupby('gene_name').agg({
        'ecc_id': 'nunique',
        'gene_type': 'first',
        'size': 'median'
    }).rename(columns={'ecc_id': 'ecc_count', 'size': 'median_size'})
    
    gene_counts = gene_counts.sort_values('ecc_count', ascending=False)
    hotspots = gene_counts.head(top_n).reset_index()
    
    print("\n" + "="*60)
    print("eccDNA HOTSPOTS (Genes with Most Overlapping eccDNA)")
    print("="*60)
    print(f"\nTop {top_n} genes:")
    print("-" * 60)
    print(f"{'Gene':<20} {'Type':<20} {'eccDNA Count':>12} {'Median Size':>12}")
    print("-" * 60)
    for _, row in hotspots.iterrows():
        print(f"{row['gene_name']:<20} {row['gene_type']:<20} {row['ecc_count']:>12,} {row['median_size']:>12,.0f}")
    
    # Visualization
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = ['steelblue' if gt == 'protein_coding' else 'coral' 
              for gt in hotspots['gene_type']]
    
    bars = ax.barh(range(len(hotspots)), hotspots['ecc_count'], color=colors)
    ax.set_yticks(range(len(hotspots)))
    ax.set_yticklabels(hotspots['gene_name'])
    ax.set_xlabel('Number of Overlapping eccDNA')
    ax.set_title(f'Top {top_n} eccDNA Hotspot Genes')
    ax.invert_yaxis()
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='steelblue', label='Protein Coding'),
                       Patch(facecolor='coral', label='Other')]
    ax.legend(handles=legend_elements, loc='lower right')
    
    # Add count labels
    for i, count in enumerate(hotspots['ecc_count']):
        ax.text(count + max(hotspots['ecc_count']) * 0.01, i, f'{count:,}', va='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'eccdna_hotspots.png')
    plt.close()
    
    return hotspots


def analyze_distance_distribution(unique_ecc: pd.DataFrame, output_dir: Path) -> dict:
    """
    Analyze distribution of distances to nearest gene.
    """
    distances = unique_ecc['distance']
    
    genic = (distances == 0).sum()
    intergenic = (distances > 0).sum()
    
    # Distance statistics for intergenic eccDNA
    intergenic_distances = distances[distances > 0]
    
    stats_dict = {
        'genic_count': genic,
        'genic_pct': genic / len(distances) * 100,
        'intergenic_count': intergenic,
        'intergenic_pct': intergenic / len(distances) * 100,
        'intergenic_median_distance': intergenic_distances.median() if len(intergenic_distances) > 0 else 0,
        'intergenic_mean_distance': intergenic_distances.mean() if len(intergenic_distances) > 0 else 0,
    }
    
    print("\n" + "="*60)
    print("DISTANCE TO NEAREST GENE ANALYSIS")
    print("="*60)
    print(f"\nGenic (distance=0):      {stats_dict['genic_count']:,} ({stats_dict['genic_pct']:.1f}%)")
    print(f"Intergenic (distance>0): {stats_dict['intergenic_count']:,} ({stats_dict['intergenic_pct']:.1f}%)")
    
    if len(intergenic_distances) > 0:
        print(f"\nIntergenic eccDNA distances:")
        print(f"  Median: {stats_dict['intergenic_median_distance']:,.0f} bp")
        print(f"  Mean:   {stats_dict['intergenic_mean_distance']:,.0f} bp")
    
    # Create visualization
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # 1. Pie chart
    ax1 = axes[0]
    ax1.pie([genic, intergenic], labels=['Genic', 'Intergenic'], 
            autopct='%1.1f%%', colors=['forestgreen', 'orange'],
            explode=(0.02, 0), startangle=90)
    ax1.set_title('Genic vs Intergenic eccDNA')
    
    # 2. Histogram of intergenic distances
    ax2 = axes[1]
    if len(intergenic_distances) > 0:
        ax2.hist(intergenic_distances, bins=50, log=True, color='orange', edgecolor='white', alpha=0.8)
        ax2.axvline(intergenic_distances.median(), color='red', linestyle='--', linewidth=2,
                   label=f'Median ({intergenic_distances.median():,.0f} bp)')
        ax2.set_xlabel('Distance to Nearest Gene (bp)')
        ax2.set_ylabel('Count (log scale)')
        ax2.set_title('Distance Distribution (Intergenic eccDNA)')
        ax2.legend()
    
    # 3. Distance categories
    ax3 = axes[2]
    if len(intergenic_distances) > 0:
        bins = [0, 1000, 5000, 10000, 50000, 100000, float('inf')]
        labels = ['<1kb', '1-5kb', '5-10kb', '10-50kb', '50-100kb', '>100kb']
        categories = pd.cut(intergenic_distances, bins=bins, labels=labels)
        cat_counts = categories.value_counts().reindex(labels)
        
        bars = ax3.bar(labels, cat_counts.values, color=plt.cm.Oranges(np.linspace(0.3, 0.9, len(labels))))
        ax3.set_xlabel('Distance Category')
        ax3.set_ylabel('Count')
        ax3.set_title('Intergenic eccDNA by Distance Category')
        ax3.tick_params(axis='x', rotation=45)
        
        # Add count labels
        for bar, count in zip(bars, cat_counts.values):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(cat_counts)*0.02, 
                    f'{count:,}', ha='center', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'distance_analysis.png')
    plt.close()
    
    return stats_dict


def size_vs_distance_analysis(unique_ecc: pd.DataFrame, output_dir: Path):
    """
    Analyze relationship between eccDNA size and distance to genes.
    """
    # Filter to intergenic for distance analysis
    intergenic = unique_ecc[unique_ecc['distance'] > 0].copy()
    
    if len(intergenic) < 10:
        print("\nInsufficient intergenic eccDNA for size vs distance analysis")
        return
    
    # Spearman correlation (robust to outliers)
    corr, p_value = stats.spearmanr(intergenic['size'], intergenic['distance'])
    
    print("\n" + "="*60)
    print("SIZE vs DISTANCE CORRELATION ANALYSIS")
    print("="*60)
    print(f"\nSpearman correlation: ρ = {corr:.4f}")
    print(f"P-value: {p_value:.2e}")
    if p_value < 0.05:
        if corr > 0:
            print("→ Significant positive correlation: larger eccDNA tend to be further from genes")
        else:
            print("→ Significant negative correlation: larger eccDNA tend to be closer to genes")
    else:
        print("→ No significant correlation between size and distance")
    
    # Create scatter plot with regression
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Subsample for visualization if too many points
    if len(intergenic) > 5000:
        sample = intergenic.sample(5000, random_state=42)
    else:
        sample = intergenic
    
    ax.scatter(sample['size'], sample['distance'], alpha=0.3, s=10, c='steelblue')
    
    # Add regression line
    z = np.polyfit(np.log10(sample['size'] + 1), np.log10(sample['distance'] + 1), 1)
    p = np.poly1d(z)
    x_line = np.logspace(np.log10(sample['size'].min()), np.log10(sample['size'].max()), 100)
    ax.plot(x_line, 10**p(np.log10(x_line)), 'r-', linewidth=2, label=f'Trend (ρ={corr:.3f})')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('eccDNA Size (bp)')
    ax.set_ylabel('Distance to Nearest Gene (bp)')
    ax.set_title('Relationship Between eccDNA Size and Distance to Genes')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'size_vs_distance.png')
    plt.close()


def generate_summary_report(stats: dict, output_dir: Path) -> None:
    """Generate a text summary report."""
    report = f"""
================================================================================
                    eccDNA ANALYSIS SUMMARY REPORT
================================================================================

Data Overview
-------------
Total unique eccDNA detected: {stats['size']['count']:,}
  - microDNA (<400bp):  {stats['size']['microDNA_count']:,} ({stats['size']['microDNA_pct']:.1f}%)
  - eccDNA (≥400bp):    {stats['size']['eccDNA_count']:,} ({100-stats['size']['microDNA_pct']:.1f}%)

Size Statistics
---------------
  Mean:   {stats['size']['mean']:,.1f} bp
  Median: {stats['size']['median']:,.1f} bp
  Std:    {stats['size']['std']:,.1f} bp
  Range:  {stats['size']['min']:,} - {stats['size']['max']:,} bp
  IQR:    {stats['size']['q25']:,.1f} - {stats['size']['q75']:,.1f} bp

Genomic Location
----------------
  Genic (overlapping genes):    {stats['distance']['genic_count']:,} ({stats['distance']['genic_pct']:.1f}%)
  Intergenic (between genes):   {stats['distance']['intergenic_count']:,} ({stats['distance']['intergenic_pct']:.1f}%)

Generated Figures
-----------------
  1. size_distribution_comprehensive.png
  2. gene_type_enrichment.png
  3. eccdna_hotspots.png
  4. distance_analysis.png
  5. size_vs_distance.png

================================================================================
"""
    
    report_path = output_dir / 'analysis_summary.txt'
    with open(report_path, 'w') as f:
        f.write(report)
    
    print(f"\n\nSummary report saved to: {report_path}")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(
        description='Extended analysis of annotated eccDNA data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('input_file', help='Path to annotated eccDNA TSV file')
    parser.add_argument('--output-dir', '-o', default='analysis_output',
                        help='Output directory for figures and reports')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print("="*60)
    print("EXTENDED eccDNA ANALYSIS")
    print("="*60)
    print(f"\nInput file: {args.input_file}")
    print(f"Output directory: {output_dir}")
    
    # Load data
    print("\nLoading data...")
    df = load_data(args.input_file)
    unique_ecc = get_unique_eccdna(df)
    print(f"Loaded {len(df):,} rows, {len(unique_ecc):,} unique eccDNA")
    
    # Run analyses
    all_stats = {}
    
    all_stats['size'] = analyze_size_distribution(unique_ecc, output_dir)
    analyze_gene_type_enrichment(df, output_dir)
    identify_hotspots(df, output_dir)
    all_stats['distance'] = analyze_distance_distribution(unique_ecc, output_dir)
    size_vs_distance_analysis(unique_ecc, output_dir)
    
    # Generate summary report
    generate_summary_report(all_stats, output_dir)
    
    print(f"\nAll figures saved to: {output_dir}/")
    print("\nAnalysis complete!")


if __name__ == '__main__':
    main()
