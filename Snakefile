rule aact_data_map:
    input:
        "data/in/combined_neuro_designs_calculated_full_v2_20240701.csv",
        "data/in/flat_ner_annotations_basic_dict_mapped_19632.csv",
        "data/in/studies_enrollment_20240717.csv"
    output:
        "data/out/all_metadata_filtered.csv",
        "data/out/general_metadata.csv",
        "data/out/trial_design_metadata.csv",
        "data/out/trial_results_reporting_metadata.csv",
        "data/out/trial_sponsorship_metadata.csv",
        "data/out/trial_enrollment_metadata.csv"
    script:
        "scripts/prepare_aact_metadata.py"

rule viz_general_metadata:
    input:
        "data/out/general_metadata.csv"
    output:
        "viz/trials_purpose_and_status.pdf",
        "viz/trials_phase_and_over_time.pdf"
    script:
        "scripts/viz_general_metadata.py"

rule viz_design_characteristics:
    input:
        "data/out/trial_design_metadata.csv",
        "data/out/trial_enrollment_metadata.csv"
    output:
        "viz/design_allocation_over_time.pdf",
        "viz/design_masking_over_time.pdf",
        "viz/design_trial_size_facilities_enrollment.pdf"
    script:
        "scripts/viz_design_characteristics.py"
