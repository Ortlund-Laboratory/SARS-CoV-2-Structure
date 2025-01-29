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
[data/PacBio_amplicons.gb]XXX<br>
GeneBank data file describing sequence features.<br>
[data/barcode_runs.csv]XXX<br>
List of Illumina barcode samples to be analyzed by the snakemake workflow.<br>
[data/processed_ccs.csv]XXX<br>
Processed PacBio CCSs, generated from our [PacBio_Library_Sequencing]XXX routine. Ensure the library is consistent with those used for the assay.<br>
[data/wildtype_sequence.fasta]XXX<br>
SARS-CoV-2 nucleocapsid wildtype sequence (Wuhan).<br>

### Sequencing Data

The workflow operates on Illumina barcode sequencing data in fastq.gz format and these files are kept compressed throughout. File location and name should match the listings given in [data/barcode_runs.csv]XXX. These files are too large to be contained in GitHub, and so are found, respectively, at:

**I SHOULD UPLOAD THESE FILES TO NCBI SRA**<br>

## Workflow

Use the `snakemake` environment:

`conda activate snakemake`

Run `snakemake` using specified number of cores:

`snakemake -j 6`

## Key Output

[results/counts/barcode_fates.csv]XXX<br>
Tally of barcodes classified and filtered according to quality.<br>
[results/counts/variant_counts.csv]XXX<br>
Tally of individual barcode counts for each sample.<br>

## Stability Score Generation & Data Visualization

Go to [scores_and_visualization]XXX for higher-/lower-stability score calculation and heatmap generation.
