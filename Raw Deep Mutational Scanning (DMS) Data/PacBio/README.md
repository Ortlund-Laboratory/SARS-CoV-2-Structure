# PacBio Library Sequencing

PacBio library sequencing is required to connect variants to their respective barcodes.

Input fastq files are too large for GitHub, and can be found at the NCBI Sequence Read Archive (SRA).

Search NCBI SRA Accession SRR28373652 or click [here](https://www.ncbi.nlm.nih.gov/sra/SRX23978660[accn]).<br>
Change line above!!! At the moment, it points to wrong SRA.

## Input Files Required

[config.yaml]XXX<br>
Configuration script controlling variables used by Jupyter notebook.<br>
[data/feature_parse_specs.yaml]XXX<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb]XXX<br>
GeneBank data file describing sequence features.<br>
[data/PacBio_runs.csv]XXX<br>
List of sequence (fastq) files to be analyzed.<br>
**results/ccs/XXX.fastq**<br>
Input circular consensus sequences (CCSs) data file.<br>
[process_ccs.ipynb]XXX<br>
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

[results/process_ccs/processed_ccs.csv]XXX<br>
List of parsed barcodes and their associated mutations.



