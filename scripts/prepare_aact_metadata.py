import pandas as pd

def main(metadata_file, ner_annotations_file, output_files):
    # Load the input files
    trial_metadata = pd.read_csv(metadata_file)
    ner_annotations = pd.read_csv(ner_annotations_file)[['nct_id']]

    # Filter and merge the dataframes
    df = trial_metadata[trial_metadata['nct_id'].isin(ner_annotations['nct_id'])]

    # Print information about the sizes of the dataframes
    print('Size trial ner_annotations: ', len(trial_metadata))
    print('Size trial metada all before filter: ', len(trial_metadata))
    print('Size trial metadata all after filter: ', len(df))
    print('Unique NCTIDs saved: ', len(set(df['nct_id'])))

    df['start_date'] = pd.to_datetime(df['start_date'])
    df['start_year'] = df['start_date'].dt.year
    df['completion_date'] = pd.to_datetime(df['completion_date'])
    df['completion_year'] = df['completion_date'].dt.year

    # Save the merged dataframe to the output file
    df.to_csv(output_files[0], index=False)

    trial_metadata = df[['nct_id','start_year', 'completion_year', 'phase', 'overall_status','primary_purpose']].drop_duplicates()
    print('Size trial_metadata: ', len(trial_metadata))
    trial_metadata.to_csv(output_files[1], index=False)

    trial_design = df[['nct_id','allocation', 'masking', 'number_of_primary_outcomes_to_measure', 'number_of_secondary_outcomes_to_measure', 'number_of_other_outcomes_to_measure','number_of_facilities', 'country', 'country_name', 'start_year']].drop_duplicates()
    print('Size trial_design: ', len(trial_design))
    trial_design.to_csv(output_files[2], index=False)

    df_results_reported = df[['nct_id', 'start_year', 'completion_year','were_results_reported', 'months_to_report_results', 'overall_status']].drop_duplicates()
    print('Size df_results_reported: ', len(df_results_reported))
    df_results_reported.to_csv(output_files[3], index=False)

    df_funding = df[['nct_id', 'start_year', 'agency_class', 'lead_or_collaborator', 'sponsor_name', 'phase']].drop_duplicates()
    print('Size df_funding: ', len(df_funding))
    df_funding.to_csv(output_files[4], index=False)

    df_country = df[['nct_id', 'country']].drop_duplicates()
    print('Size df_country: ', len(df_country))
    df_country.to_csv(output_files[5], index=False)


if __name__ == "__main__":
    if "snakemake" in globals():
        metadata_file = snakemake.input[0]
        ner_annotations_file = snakemake.input[1]
        output_files = snakemake.output
        main(metadata_file, ner_annotations_file, output_files)
    else:
        raise RuntimeError("This script should be run using Snakemake.")
