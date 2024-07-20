import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def viz_results_reporting(trial_metadata, output_file):

    df_results_reported = trial_metadata[trial_metadata['overall_status']=='Completed']
    df_results_reported = df_results_reported[df_results_reported['completion_year']<2022]

    year_to_use = 'completion_year'

    # Group by completion_year and were_results_reported, then count unique nct_ids
    results_over_time = df_results_reported.groupby([year_to_use, 'were_results_reported'])['nct_id'].nunique().unstack(fill_value=0)

    # Calculate the proportion of each type
    results_proportion_over_time = results_over_time.div(results_over_time.sum(axis=1), axis=0)

    # Filter data for reported results
    reported_results = df_results_reported[df_results_reported['were_results_reported'] == True]

    # Calculate average months to report results and standard deviation
    average_months_to_report = reported_results.groupby(year_to_use)['months_to_report_results'].mean()
    std_months_to_report = reported_results.groupby(year_to_use)['months_to_report_results'].std()

    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

    # Subplot 1: Stacked bar chart for results reported/results pending
    bottom = np.zeros(len(results_proportion_over_time))
    colors = {False: '#DDCC77', True: '#44AA99'}
    for reported in [False, True]:
        label = 'Results Reported' if reported else 'Results Pending'
        bars = ax1.bar(results_proportion_over_time.index, 
                    results_proportion_over_time[reported], 
                    bottom=bottom, 
                    label=label, 
                    color=colors[reported], zorder=2)
        bottom += results_proportion_over_time[reported]

        # Add labels to each segment
        for bar in bars:
            height = bar.get_height()
            if height > 0.05:
                ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_y() + height / 2, f'{height*100:.0f}%', 
                        ha='center', va='center', fontsize=10, rotation=90)

    ax1.set_xlabel("Trial " + year_to_use.replace("_", " ").capitalize(), fontsize=14)
    #ax1.set_ylabel('Proportion of Trials', fontsize=14)
    ax1.set_title('Proportion Reported Results for Trials Completed in that Year', fontsize=15)
    ax1.legend(loc='lower left', title='', fontsize=13)
    ax1.grid(linestyle='--', alpha=0.6, zorder=1)
    ax1.tick_params(axis='x', labelsize=14)
    ax1.tick_params(axis='y', labelsize=14)
    ax1.set_ylim(0, 1)
    ax1.set_yticklabels(['{:.0f}%'.format(x * 100) for x in ax1.get_yticks()])
    ax1.text(-0.03, 1.08, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')


    # Subplot 2: Line chart for months to report results
    ax2.plot(average_months_to_report.index, average_months_to_report, marker='o', label='Avg Months to Report Results', color='#882255')
    ax2.fill_between(average_months_to_report.index, average_months_to_report - std_months_to_report, average_months_to_report + std_months_to_report, color='lightgrey', alpha=0.3)
    ax2.set_xlabel("Trial " + year_to_use.replace("_", " ").capitalize(), fontsize=14)
    #ax2.set_ylabel('Average Months to Report Results', fontsize=14)
    ax2.set_title('Average Number of Months to Report Results for Trials Completed in that Year', fontsize=15)
    #ax2.legend(loc='upper right', title='', fontsize=13)
    ax2.grid(linestyle='--', alpha=0.6, zorder=1)
    ax2.tick_params(axis='x', labelsize=14)
    ax2.tick_params(axis='y', labelsize=14)
    ax2.text(-0.03, 1.08, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    plt.tight_layout()

    # Save the figure as a PDF
    plt.savefig(output_file)

def main(metadata_file, output_file):
    # Load the input files
    trial_metadata = pd.read_csv(metadata_file)
    viz_results_reporting(trial_metadata, output_file)

if __name__ == "__main__":
    if "snakemake" in globals():
        metadata_file = snakemake.input[0]
        output_file = snakemake.output[0]
        main(metadata_file, output_file)
    else:
        raise RuntimeError("This script should be run using Snakemake.")
