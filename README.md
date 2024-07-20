# Triomics_Pipeline

This project is a data analysis pipeline designed to process and visualize clinical trials data using Snakemake. 
The pipeline reads input data files, processes them to extract relevant metadata, and generates visualizations to help understand various aspects of clinical trials. 
The workflow is automated using Snakemake to ensure reproducibility and efficiency.

Snakemake is a workflow management system that enables the creation of reproducible and scalable data analyses. 
For more information, visit the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

## Setup
To set up the project, follow these steps:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/Ineichen-Group/Triomics_Pipeline.git
   cd Triomics-Pipeline
   ```
2. **Create a conda environment:**
   1. Option 1 - from the yml file:
    ```bash
    conda env create -f environment.yml
    conda activate neurotrial-analysis
   ```
   2. Clean environment in which later you install the dependencies with pip
   ```bash
    conda create --name clinical-trials-analysis python=3.8
    conda activate clinical-trials-analysis
   ```
   ```bash
    pip install -r requirements.txt
   ```
3. **Ensure all necessary directories are created:**
    ```bash
    mkdir -p data/in data/out viz
   ```

## Workflow Description

The pipeline consists of several rules that define the workflow. Each rule specifies input files, output files, and a script to process the data. Below is a description of each rule in the pipeline:

### Rule: aact_data_map

**Description:** This rule processes the input data files to generate various metadata files.
