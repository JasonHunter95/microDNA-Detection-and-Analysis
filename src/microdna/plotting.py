#!/usr/bin/env python3
"""
plotting.py - Shared plotting utilities for eccDNA analysis.

Provides reusable visualization functions for eccDNA length and score
distributions with consistent styling.
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_distribution(
    data: pd.DataFrame,
    column: str,
    title: str,
    xlabel: str,
    ylabel: str = "Frequency",
    log_scale: bool = True,
    figsize: tuple[int, int] = (12, 6),
) -> plt.Figure:
    """
    Plot the distribution of a column with KDE overlay and statistical markers.

    Args:
        data: DataFrame containing the data.
        column: Column name to plot.
        title: Title of the plot.
        xlabel: Label for the x-axis.
        ylabel: Label for the y-axis.
        log_scale: Whether to use log scale on x-axis.
        figsize: Figure size as (width, height).

    Returns:
        The matplotlib Figure object.
    """
    fig = plt.figure(figsize=figsize)

    # Filter out zero or negative values before log transform
    data_to_plot = data[data[column] > 0]

    sns.histplot(
        data=data_to_plot,
        x=column,
        bins="auto",
        kde=True,
        log_scale=log_scale,
    )

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Add statistical markers
    plt.axvline(
        data[column].mean(), color="r", linestyle="--", label="Mean"
    )
    plt.axvline(
        data[column].median(), color="g", linestyle="--", label="Median"
    )
    plt.axvline(
        data[column].quantile(0.75), color="b", linestyle="--", label="75th Percentile"
    )
    plt.axvline(
        data[column].quantile(0.25), color="y", linestyle="--", label="25th Percentile"
    )
    plt.legend()

    return fig
