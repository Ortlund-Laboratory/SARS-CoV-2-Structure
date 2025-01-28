# PacBio Library Sequencing

PacBio library sequencing is required to connect variants to their respective barcodes.

Input fastq files are too large for GitHub, and can be found at the NCBI Sequence Read Archive (SRA).

Search NCBI SRA Accession SRR28373652 or click [here](https://www.ncbi.nlm.nih.gov/sra/SRX23978660[accn]).<br>
Change line above!!! At the moment, it points to wrong SRA.

## Input Files Required

[config.yaml](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/PacBio/config.yaml)<br>
Configuration script controlling variables used by Jupyter notebook.<br>
[data/feature_parse_specs.yaml](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/PacBio/data/feature_parse_specs.yaml)<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/PacBio/data/PacBio_amplicons.gb)<br>
GeneBank data file describing sequence features.<br>
[data/PacBio_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/PacBio/data/PacBio_runs.csv)<br>
List of sequence (fastq) files to be analyzed.<br>
**results/ccs/XXX.fastq**<br>
Input circular consensus sequences (CCSs) data file.<br>
[process_ccs.ipynb](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/PacBio/process_ccs.ipynb)<br>
Jupyter notebook for extracting barcodes from CCSs and matching to variants.<br>

## Setup

Ensure the fastq file(s) listed in PacBio_runs.csv are present in results/ccs BEFORE executing the workflow.

## Workflow

Use the `snakemake` environment:

`conda activate snakemake`

Run `jupyter`:

`jupyter notebook process_ccs.ipynb`

**NOTE:** Some cells of the `jupyter` notebook may fail to execute, in the absence of ccs summary files. These are not required, and so the user can skip to the next cell.

## Key Output

[results/process_ccs/processed_ccs.csv.gz]XXX<br>
List of parsed barcodes and their associated mutations.



