import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def viz_purpose_status(trial_metadata, output_file):
    # Keeping only unique nct_id, primary_purpose pairs
    unique_pairs_purpose = trial_metadata[['nct_id', 'primary_purpose']].drop_duplicates()

    # Counting the number of nct_ids per primary purpose type
    purpose_type_counts = unique_pairs_purpose['primary_purpose'].value_counts()
    purpose_type_counts = purpose_type_counts.sort_values(ascending=True)

    # Keeping only unique nct_id, overall_status pairs
    unique_pairs_status = trial_metadata[['nct_id', 'overall_status']].drop_duplicates()

    # Counting the number of nct_ids per overall status type
    status_type_counts = unique_pairs_status['overall_status'].value_counts()
    status_type_counts = status_type_counts.sort_values(ascending=True)

    # Create a figure with two horizontal bar charts
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))

    # Plot for Primary Purpose
    ax1.barh(purpose_type_counts.index, purpose_type_counts, color='lightgrey', zorder=2)
    for i, v in enumerate(purpose_type_counts):
        ax1.text(v + 50, i, str(v), va='center', color='black', fontsize=12)
    ax1.set_title('Primary Trial Purpose', fontsize=14)
    ax1.set_xlabel('Count of Unique Trials', fontsize=14)
    ax1.set_xlim(0, max(purpose_type_counts)+1500)  # Adjusting the x limits for visibility
    ax1.grid(axis='x', linestyle='--', alpha=0.6, zorder=1)
    ax1.set_xlim(0, max(purpose_type_counts)+1800) # Adjusted to max count for relevancy
    ax1.tick_params(axis='both', labelsize=13)  # Increase tick label size
    ax1.text(-0.03, 1.05, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Plot for Overall Status
    ax2.barh(status_type_counts.index, status_type_counts, color='lightgrey', zorder=2)
    for i, v in enumerate(status_type_counts):
        ax2.text(v + 50, i, str(v), va='center', color='black', fontsize=12)
    ax2.set_title('Overall Trial Status', fontsize=14)
    ax2.set_xlabel('Count of Unique Trials', fontsize=14)
    ax2.set_xlim(0, max(status_type_counts)+1000)  # Adjusting the x limits for visibility
    ax2.grid(axis='x', linestyle='--', alpha=0.6, zorder=1)
    ax2.set_xlim(0, max(status_type_counts)+1500) # Adjusted to max count for relevancy
    ax2.tick_params(axis='both', labelsize=13)  # Increase tick label size
    ax2.text(-0.03, 1.05, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    # Adjust layout and display the plots
    plt.tight_layout()

    # Optionally save the figure to a local folder
    fig.savefig(output_file)

def viz_phase_year_growth(trial_metadata, output_file):
    # Define phase order and filter/preprocess the data
    phase_order = [
        'Early Phase 1',
        'Phase 1',
        'Phase 1/2',
        'Phase 2',
        'Phase 2/3',
        'Phase 3',
        'Phase 4',
        'Not Applicable'
    ]

    # Define custom colors for each phase
    phase_colors = {
        'Early Phase 1': '#882255',
        'Phase 1': '#AA4499',
        'Phase 1/2': '#CC6677',
        'Phase 2': '#DDCC77',
        'Phase 2/3': '#88CCEE',
        'Phase 3': '#44AA99',
        'Phase 4': '#117733',
        'Not Applicable': '#332288'
    }

    # Filter and count phases and statuses
    unique_pairs_phase = trial_metadata[['nct_id', 'phase', 'overall_status']].drop_duplicates()
    unique_pairs_phase['phase'] = unique_pairs_phase['phase'].str.replace('/Phase ', '/', regex=False)

    phase_type_counts = unique_pairs_phase['phase'].value_counts().reindex(phase_order, fill_value=0)
    completed_count = unique_pairs_phase[unique_pairs_phase['overall_status'] == 'Completed']['phase'].value_counts().reindex(phase_order, fill_value=0)
    completed_proportion = (completed_count / phase_type_counts * 100).fillna(0)

    # Prepare data for the time series plot
    filtered_data = trial_metadata[['nct_id', 'phase', 'start_year']][trial_metadata['start_year'] < 2024].drop_duplicates()
    trial_counts = filtered_data.groupby(['phase', 'start_year']).size().unstack(fill_value=0)
    total_trials_per_year = trial_counts.sum(axis=0)

    # Create subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 5))

    # First subplot - Completion by phase
    total_bars = ax1.bar(phase_type_counts.index, phase_type_counts, color='lightgrey', zorder=2, label='Total Trials')
    completed_bars = ax1.bar(completed_count.index, completed_count, color='darkgrey', zorder=2, label='Completed Trials')

    # Label the bars
    for bar, prop in zip(completed_bars, completed_proportion):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height(), f'{prop:.1f}%', ha='center', va='bottom', fontsize=10)

    for bar in total_bars:
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height(), f'{int(bar.get_height())}', ha='center', va='bottom', fontsize=12)

    ax1.set_xlabel('Trial Phase',fontsize=14)
    ax1.set_title('Trial Completion by Phase',fontsize=14)
    ax1.legend()
    ax1.grid(axis='y', linestyle='--', alpha=0.6)
    ax1.text(-0.03, 1.05, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    ax1.tick_params(axis='both', labelsize=12)  # Increase tick label size

    # Second subplot - Trials over time
    ax2.fill_between(total_trials_per_year.index, total_trials_per_year, color='lightgray', alpha=0.3, label='Total Trials')
    for phase in phase_order:
        if phase in trial_counts.index:
            ax2.plot(trial_counts.columns, trial_counts.loc[phase], label=phase, color=phase_colors[phase])

    ax2.set_xlabel('Start Year', fontsize=14)
    ax2.set_title('Count of Trials Started by Phase and Year',fontsize=14)
    ax2.legend(loc='upper left')
    ax2.grid(linestyle='--', alpha=0.6, zorder=1)
    ax2.set_xticks(np.arange(min(trial_counts.columns), max(trial_counts.columns)+1, 5))
    ax2.text(-0.03, 1.05, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    ax2.tick_params(axis='both', labelsize=13)  # Increase tick label size
    # Adjust layout and display the plots
    plt.tight_layout()

    # Optionally save the figure to a local folder
    fig.savefig(output_file)


def main(metadata_file, output_files):
    # Load the input files
    trial_metadata = pd.read_csv(metadata_file)
    viz_purpose_status(trial_metadata, output_files[0])
    viz_phase_year_growth(trial_metadata, output_files[1])


if __name__ == "__main__":
    if "snakemake" in globals():
        metadata_file = snakemake.input[0]
        output_files = snakemake.output
        main(metadata_file, output_files)
    else:
        raise RuntimeError("This script should be run using Snakemake.")
