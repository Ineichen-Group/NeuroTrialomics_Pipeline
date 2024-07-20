import pandas as pd

def classify_enrollment(enrollment):
    if 0 <= enrollment <= 10:
        return "0–10"
    elif 11 <= enrollment <= 50:
        return "11–50"
    elif 51 <= enrollment <= 100:
        return "51–100"
    elif 101 <= enrollment <= 1000:
        return "101–1,000"
    elif enrollment > 1000:
        return ">1,000"
    else:
        return "not reported"
    
def main(metadata_file, ner_annotations_file, enrollment_file, output_files):
    # Load the input files
    trial_metadata = pd.read_csv(metadata_file)
    ner_annotations = pd.read_csv(ner_annotations_file)[['nct_id']]
    aact_baseline_counts = pd.read_csv(enrollment_file)

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
    # Fill NA values in all columns except 'start_year'
    columns_to_fill = trial_design.columns.difference(['start_year'])
    trial_design[columns_to_fill] = trial_design[columns_to_fill].fillna('not reported')
            
    #trial_design.fillna('not reported', inplace=True)
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

    trials_ids_unique = df[['nct_id','start_year']].drop_duplicates()
    trials_with_participants = trials_ids_unique.merge(aact_baseline_counts, how='left', on='nct_id')
    trials_with_participants['enrollment_class'] = trials_with_participants['enrollment'].apply(classify_enrollment)
    print('Size trials_with_participants: ', len(trials_with_participants))
    trials_with_participants.to_csv(output_files[6], index=False)


if __name__ == "__main__":
    if "snakemake" in globals():
        metadata_file = snakemake.input[0]
        ner_annotations_file = snakemake.input[1]
        enrollment_file = snakemake.input[2]
        output_files = snakemake.output
        main(metadata_file, ner_annotations_file, enrollment_file, output_files)
    else:
        raise RuntimeError("This script should be run using Snakemake.")
