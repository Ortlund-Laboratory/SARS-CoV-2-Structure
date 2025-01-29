# Workflow to Calculate Stability Scores

SARS-CoV-2 nucleocapsid variant stability according to Myc selection was investigated, and reference, high-Myc and low-Myc populations collected for sequencing. These sequences can be processed to give scores for both more and less stable variants.

## Input Files Required

[SnakeFile](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/Snakefile)<br>
Gives overall instructions for the `snakemake` workflow.<br>
[config.yaml](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/config.yaml)<br>
Configuration script controlling variables used by Jupyter notebooks.<br>
[build_variants.ipynb](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/build_variants.ipynb)<br>
Builds a barcode variant table based on the data from the processed PacBio CCSs.<br>
[R2_to_R1.py](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/R2_to_R1.py)<br>
Converts barcodes located at the R2 end to the R1 end by taking the reverse complement. This allows the barcodes to be read and parsed correctly by the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) algorithm.<br>
[count_variants.ipynb](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/count_variants.ipynb)<br>
Counts the number of times a barcode (and by extension a variant) appears in each Illumina barcode sequencing sample.<br>
[scripts/run_nb.py](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scripts/run_nb.py)<br>
Runs Jupyter notebooks and creates Markdown output.<br>
[data/feature_parse_specs.yaml](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/data/feature_parse_specs.yaml)<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/data/PacBio_amplicons.gb)<br>
GeneBank data file describing sequence features.<br>
[data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/data/barcode_runs.csv)<br>
List of Illumina barcode samples to be analyzed by the snakemake workflow.<br>
[data/processed_ccs.csv.gz](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/data/processed_ccs.csv.gz)<br>
Processed PacBio CCSs, generated from our [PacBio](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/tree/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/PacBio) routine. Ensure the library is consistent with those used for the assay.<br>
[data/wildtype_sequence.fasta](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/data/wildtype_sequence.fasta)<br>
SARS-CoV-2 nucleocapsid wildtype sequence (Wuhan).<br>

### Sequencing Data

The workflow operates on Illumina barcode sequencing data in fastq.gz format and these files are kept compressed throughout. File location and name should match the listings given in [data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/data/barcode_runs.csv). These files are too large to be contained in GitHub, and so are found, respectively, at:

**p22064-s003_220225_03_S150_L002_R2_001.fastq.gz**<br>
Search NCBI BioProject PRJNA1216977 (BioSample SAMN46474440) or click [here](https://www.ncbi.nlm.nih.gov/bioproject/1216977).

**p22064-s004_220225_04_S151_L002_R2_001.fastq.gz**<br>
Search NCBI BioProject PRJNA1216977 (BioSample SAMN46474441) or click [here](https://www.ncbi.nlm.nih.gov/bioproject/1216977).

**p22064-s005_220225_05_S152_L002_R2_001.fastq.gz**<br>
Search NCBI BioProject PRJNA1216977 (BioSample SAMN46474442) or click [here](https://www.ncbi.nlm.nih.gov/bioproject/1216977).

## Workflow

Use the `snakemake` environment:

`conda activate snakemake`

Run `snakemake` using specified number of cores:

`snakemake -j 6`

## Key Output

[results/counts/barcode_fates.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/results/counts/barcode_fates.csv)<br>
Tally of barcodes classified and filtered according to quality.<br>
[results/counts/variant_counts.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/results/counts/variant_counts.csv)<br>
Tally of individual barcode counts for each sample.<br>

## Stability Score Generation & Data Visualization

Go to [scores_and_visualization]XXX for higher-/lower-stability score calculation and heatmap generation.
