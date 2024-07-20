# Triomics_Pipeline

This project is a data analysis pipeline designed to process and visualize clinical trials data using Snakemake. The pipeline reads input data files, processes them to extract relevant metadata, and generates visualizations to help understand various aspects of clinical trials. The workflow is automated using Snakemake to ensure reproducibility and efficiency.

## Workflow Description

The pipeline consists of several rules that define the workflow. Each rule specifies input files, output files, and a script to process the data. Below is a description of each rule in the pipeline:

### Rule: aact_data_map

**Description:** This rule processes the input data files to generate various metadata files.
