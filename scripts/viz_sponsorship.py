import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import contextily as ctx
import geopandas as gpd
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import ScalarMappable

def viz_sponsorship_lead(trial_metadata, output_file):
    df_filtered = trial_metadata[trial_metadata['lead_or_collaborator']=='lead']
    df_filtered['start_year'] = df_filtered['start_year'].astype(int)

    # Define color map
    labels = ['UNKNOWN', 'INDIV', 'NETWORK', 'OTHER_GOV', 'FED', 'NIH', 'HOSPITAL', 'OTHER', 'UNIVERSITY', 'INDUSTRY']
    colors = ['#d3d3d3', '#000000', '#999999', '#56B4E9', '#D55E00', '#F0E442', '#CC79A7', '#009E73', '#0072B2', '#E69F00']
    color_map = dict(zip(labels, colors))

    default_color = '#FFFFFF'  # White, change as needed

    # Creating a pivot table with counts per agency_class and year
    pivot_table_counts = df_filtered.pivot_table(index='start_year', columns='agency_class', values='nct_id', aggfunc='count', fill_value=0)
    pivot_table_percentage = pivot_table_counts.divide(pivot_table_counts.sum(axis=1), axis=0) * 100
    overall_distribution = pivot_table_counts.sum(axis=0).sort_values(ascending=True)

    # Reorder pivot table columns to match the order of overall_distribution
    pivot_table_counts = pivot_table_counts[overall_distribution.index]
    pivot_table_percentage = pivot_table_percentage[overall_distribution.index]

    # Plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

    # Plot 1 - Overall Distribution (now plot A)
    overall_distribution.plot(kind='barh', ax=ax1, color=[color_map.get(x, default_color) for x in overall_distribution.index], zorder=2)
    ax1.set_title('Overall Distribution of Lead Funding Sources since year 2000', fontsize=14)
    ax1.set_xlabel('Count of Unique Trials', fontsize=14)
    ax1.text(-0.03, 1.05, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    ax1.set_ylabel('')
    ax1.grid(axis='x', linestyle='--', alpha=0.6, zorder=1)

    # Adding text labels to the bars
    for index, value in enumerate(overall_distribution):
        ax1.text(value, index, f' {int(value)}', va='center', ha='left')
    ax1.set_xlim(0, max(overall_distribution) + 500)
    ax1.tick_params(axis='both', labelsize=13)  # Increase tick label size

    # Plot 2 - Yearly Distribution (now plot B)
    pivot_table_percentage.plot(kind='bar', stacked=True, ax=ax2, color=[color_map.get(x, default_color) for x in pivot_table_counts.columns])
    ax2.set_title('Percentage of Total Trials by Lead Funding Source since year 2000', fontsize=14)
    ax2.set_xlabel('Trial Start Year', fontsize=14)
    ax2.legend().set_visible(False)
    ax2.tick_params(axis='both', labelsize=13)  # Increase tick label size

    ax2.text(-0.03, 1.05, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    # Adding text inside the bars
    for bars_stack in ax2.containers:
        ax2.bar_label(bars_stack, labels=[f'{v:.0f}%' if v > 10 else '' for v in bars_stack.datavalues], label_type='center', fontsize=8)

    # Synchronize legend for both plots
    handles, labels = ax2.get_legend_handles_labels()
    ax1.legend(handles[::-1], labels[::-1], title='Funding Source Class', loc='lower right')

    plt.tight_layout()

    # Save the figure as a PDF
    plt.savefig(output_file)

def main(metadata_file, output_file):
    # Load the input files
    trial_metadata = pd.read_csv(metadata_file)
    viz_sponsorship_lead(trial_metadata, output_file)

if __name__ == "__main__":
    if "snakemake" in globals():
        metadata_file = snakemake.input[0]
        output_file = snakemake.output[0]
        main(metadata_file, output_file)
    else:
        raise RuntimeError("This script should be run using Snakemake.")
