import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def viz_allocation(trial_design, output_file):
    unique_rows = trial_design[['nct_id', 'allocation']].drop_duplicates()

    # Counting the number of nct_ids per phase type
    allocation_counts = unique_rows['allocation'].value_counts()
    allocation_counts = allocation_counts.sort_values(ascending=True)

    # Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))

    # Plot for Allocation (Plot A)
    bars_0 = axs[0].barh(allocation_counts.index, allocation_counts, color='lightgrey', zorder=2)
    for bar in bars_0:
        width = bar.get_width()
        axs[0].text(width, bar.get_y() + bar.get_height() / 2, f'{width}', va='center', fontsize=12)
    axs[0].tick_params(axis='y', labelsize=14)
    axs[0].tick_params(axis='x', labelsize=14)
    axs[0].grid(axis='x', linestyle='--', alpha=0.6, zorder=1)
    axs[0].set_xlabel('Count of Unique Trials', fontsize=14)
    axs[0].set_title('Allocation', fontsize=14)
    axs[0].set_xlim(0, max(allocation_counts) + 1500)  # Adjusted to max count for relevancy
    axs[0].text(-0.02, 1.09, 'A', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    ### PLOT B
    allocation_colors = {'Randomized': '#44AA99',
    'Non-Randomized': '#CC6677',
    'not reported': '#DDCC77'}
    allocation_types = ['Non-Randomized', 'not reported', 'Randomized']

    # Keeping only unique nct_id, allocation pairs
    unique_pairs_allocation = trial_design[['nct_id', 'allocation', 'start_year']].drop_duplicates()

    # Group by start_year and allocation, then count unique nct_ids
    allocation_over_time = unique_pairs_allocation.groupby(['start_year', 'allocation']).size().unstack(fill_value=0)

    # Calculate the proportion of each allocation type from all trials for each year
    allocation_proportion_over_time = allocation_over_time.div(allocation_over_time.sum(axis=1), axis=0)

    # Plot for Allocation Over Time (Plot B)
    bottom = np.zeros(len(allocation_proportion_over_time))
    for allocation_type in allocation_types:
        bars = axs[1].bar(allocation_proportion_over_time.index, 
                        allocation_proportion_over_time[allocation_type], 
                        bottom=bottom, 
                        label=allocation_type,
                        color=allocation_colors.get(allocation_type, 'gray'), zorder=2)  # Use custom color or gray if not specified
        bottom += allocation_proportion_over_time[allocation_type]

        # Add labels to each segment
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                axs[1].text(bar.get_x() + bar.get_width() / 2, bar.get_y() + height / 2, f'{height:.0%}', ha='center', va='center', fontsize=10, rotation=90)

    axs[1].set_xlabel('Trial Start Year', fontsize=14)
    axs[1].set_title('Proportion of Reported Allocation Over Time', fontsize=14)
    # Set y-axis limits and labels
    axs[1].set_ylim(0, 1)
    axs[1].set_yticklabels(['{:.0f}%'.format(x * 100) for x in axs[1].get_yticks()])

    # Sorting the legend handles
    handles, labels = axs[1].get_legend_handles_labels()
    allocation_type_totals = allocation_over_time.sum().to_dict()
    sorted_handles_labels = sorted(zip(handles, labels), key=lambda x: allocation_type_totals[x[1]], reverse=True)
    sorted_handles, sorted_labels = zip(*sorted_handles_labels)
    sorted_labels = [label.capitalize() for label in sorted_labels]
    axs[1].legend(sorted_handles, sorted_labels, fontsize=12, loc='upper left')

    axs[1].grid(linestyle='--', alpha=0.6, zorder=1)
    axs[1].tick_params(axis='x', labelsize=14)
    axs[1].tick_params(axis='y', labelsize=13)
    axs[1].text(-0.02, 1.09, 'B', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Adjust layout
    plt.tight_layout()

    # Save the plot to a local folder
    plt.savefig(output_file)


def viz_masking(trial_design, output_file):
    # Prepare masking data
    unique_pairs_masking = trial_design[['nct_id', 'masking', 'start_year']].drop_duplicates()
    masking_counts = unique_pairs_masking['masking'].value_counts().sort_values(ascending=True)

    # Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(15, 5))

    # Plot for Masking Frequency (Plot A)
    bars_0 = axs[0].barh(masking_counts.index, masking_counts, color='lightgrey', zorder=2)
    for bar in bars_0:
        width = bar.get_width()
        axs[0].text(width, bar.get_y() + bar.get_height() / 2, f'{width}', va='center', fontsize=12)
    axs[0].tick_params(axis='y', labelsize=14)
    axs[0].tick_params(axis='x', labelsize=14)
    axs[0].grid(axis='x', linestyle='--', alpha=0.6, zorder=1)
    axs[0].set_xlabel('Count of Unique Trials', fontsize=14)
    axs[0].set_title('Masking', fontsize=14)
    axs[0].set_xlim(0, max(masking_counts) + 630)  # Adjusted to max count for relevancy
    axs[0].text(-0.02, 1.08, 'A', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Plot for Change in Masking Over Time (Plot B)

    # Prepare masking over time data
    masking_over_time = unique_pairs_masking.groupby(['start_year', 'masking']).size().unstack(fill_value=0)

    # Calculate the proportion of each masking type from all trials for each year
    masking_proportion_over_time = masking_over_time.div(masking_over_time.sum(axis=1), axis=0)
    masking_types = ['not reported', 'None (Open Label)', 'Single', 'Double', 'Triple', 'Quadruple']

    masking_colors = {'None (Open Label)': '#88CCEE',
    'Quadruple': '#CC6677',
    'not reported': '#DDCC77',
    'Double': '#117733',
    'Triple': '#AA4499',
    'Single': '#44AA99'}

    allocation_types = ['Non-Randomized', 'not reported', 'Randomized']

    bottom = np.zeros(len(masking_proportion_over_time))
    for masking_type in masking_types:
        bars = axs[1].bar(masking_proportion_over_time.index, 
                        masking_proportion_over_time[masking_type], 
                        bottom=bottom, 
                        label=masking_type, color=masking_colors.get(masking_type, 'gray'), zorder=2)
        bottom += masking_proportion_over_time[masking_type]

        # Add labels to each segment
        for bar in bars:
            height = bar.get_height()
            #print(height)
            if height > 0.055:
                axs[1].text(bar.get_x() + bar.get_width() / 2, bar.get_y() + height / 2, f'{height:.0%}', 
                            ha='center', va='center', fontsize=10, rotation=90) 

    axs[1].set_xlabel('Trial Start Year', fontsize=14)
    axs[1].set_title('Proportion of Reported Masking Over Time', fontsize=14)

    # Set y-axis limits and labels
    axs[1].set_ylim(0, 1)
    axs[1].set_yticklabels(['{:.0f}%'.format(x * 100) for x in axs[1].get_yticks()])

    # Sort and place the legend outside the plot area
    handles, labels = axs[1].get_legend_handles_labels()
    sorted_handles_labels = sorted(zip(handles, labels), key=lambda x: masking_types[::-1].index(x[1]))
    sorted_handles, sorted_labels = zip(*sorted_handles_labels)
    axs[1].legend(sorted_handles, [label.capitalize() for label in sorted_labels], fontsize=12, loc='upper left', bbox_to_anchor=(1.00, 1))

    axs[1].grid(linestyle='--', alpha=0.6, zorder=1)
    axs[1].tick_params(axis='x', labelsize=14)
    axs[1].tick_params(axis='y', labelsize=13)
    axs[1].text(-0.02, 1.08, 'B', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Adjust layout
    plt.tight_layout()

    # Save the plot to a local folder
    plt.savefig(output_file)

def viz_facilities_enrollment(trial_design, trials_with_participants, output_file):
    labels = ['1', '2-10', '11-20', '21-30', '31-40', '41-50', '>50']

    # Prepare masking data
    unique_pairs_facilities = trial_design[['nct_id', 'binned_facilities', 'number_of_facilities']].drop_duplicates()

    # Count occurrences in each bin including NaN for 'not reported'
    bin_counts_facilities = unique_pairs_facilities['binned_facilities'].value_counts().reindex(labels + [np.nan]).fillna(0)
    bin_counts_facilities[np.nan] = len(unique_pairs_facilities[unique_pairs_facilities['number_of_facilities'].isna()])

    # Count the occurrences of each enrollment class
    enrollment_class_counts = trials_with_participants['enrollment_class'].value_counts().sort_index()

    # Define the order of the categories
    category_order_enrollment = ["0–10", "11–50", "51–100", "101–1,000", ">1,000"]

    # Reorder the counts according to the specified category order
    enrollment_class_counts = enrollment_class_counts.reindex(category_order_enrollment, fill_value=0)

    # Create subplots
    fig, axs = plt.subplots(1, 2, figsize=(16, 5))

    # Plot for number_of_facilities
    bin_counts_facilities = bin_counts_facilities.rename(index={np.nan: 'N.A.'})
    bars_facilities = axs[0].bar(bin_counts_facilities.index.astype(str), bin_counts_facilities, color='lightgrey', zorder=3)
    for i, value in enumerate(bin_counts_facilities):
        axs[0].text(i, value + 0.1, str(int(value)), ha='center', va='bottom')

    axs[0].tick_params(axis='y', labelsize=14)
    axs[0].tick_params(axis='x', labelsize=14)
    axs[0].grid(axis='y', linestyle='--', alpha=0.6, zorder=1)
    axs[0].set_xlabel('Number of Facilities', fontsize=14)
    axs[0].set_title('Number of Facilities', fontsize=14)
    #plt.xticks(rotation=45)
    axs[0].text(-0.01, 1.08, 'A', transform=axs[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Plot for enrollment
    bars_enrollment = axs[1].bar(enrollment_class_counts.index, enrollment_class_counts.values, color='lightgrey', zorder=3)
    for i, value in enumerate(enrollment_class_counts):
        axs[1].text(i, value + 0.1, str(int(value)), ha='center', va='bottom')

    axs[1].tick_params(axis='y', labelsize=14)
    axs[1].tick_params(axis='x', labelsize=14)
    axs[1].grid(axis='y', linestyle='--', alpha=0.6, zorder=1)
    axs[1].set_xlabel('Number of Participants (Actual or Anticipated)', fontsize=14)
    axs[1].set_title('Number of Enrolled Participants', fontsize=14)
    axs[1].text(-0.01, 1.08, 'B', transform=axs[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    axs[1].set_ylim(0, max(enrollment_class_counts) + 1000) # Adjusted to max count for relevancy

    plt.tight_layout()

    plt.savefig(output_file)


def main(metadata_file, enrollment_file, output_files):
    # Load the input files
    trial_metadata = pd.read_csv(metadata_file)
    trial_enrollment = pd.read_csv(enrollment_file)
    viz_allocation(trial_metadata, output_files[0])
    viz_masking(trial_metadata, output_files[1])
    viz_facilities_enrollment(trial_metadata, trial_enrollment, output_files[2])


if __name__ == "__main__":
    if "snakemake" in globals():
        metadata_file = snakemake.input[0]
        enrollment_file = snakemake.input[1]
        output_files = snakemake.output
        main(metadata_file, enrollment_file, output_files)
    else:
        raise RuntimeError("This script should be run using Snakemake.")
