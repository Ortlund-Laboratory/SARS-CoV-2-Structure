# Workflow to Calculate Initial Expression of Variants

SARS-CoV-2 nucleocapsid variant stability according to initial GFP<sup>+</sup>Myc<sup>+</sup> selection was investigated, with GFP<sup>+</sup> and GFP<sup>+</sup>Myc<sup>+</sup> populations collected for sequencing. These sequences can be processed to give scores indicating level of expression.

## Input Files Required

[SnakeFile](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/Snakefile)<br>
Gives overall instructions for the `snakemake` workflow.<br>
[config.yaml](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/config.yaml)<br>
Configuration script controlling variables used by Jupyter notebooks.<br>
[build_variants.ipynb](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/build_variants.ipynb)<br>
Builds a barcode variant table based on the data from the processed PacBio CCSs.<br>
[R2_to_R1.py](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/R2_to_R1.py)<br>
Converts barcodes located at the R2 end to the R1 end by taking the reverse complement. This allows the barcodes to be read and parsed correctly by the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) algorithm.<br>
[count_variants.ipynb](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/count_variants.ipynb)<br>
Counts the number of times a barcode (and by extension a variant) appears in each Illumina barcode sequencing sample.<br>
[scripts/run_nb.py](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/scripts/run_nb.py)<br>
Runs Jupyter notebooks and creates Markdown output.<br>
[data/feature_parse_specs.yaml](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/data/feature_parse_specs.yaml)<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/data/PacBio_amplicons.gb)<br>
GeneBank data file describing sequence features.<br>
[data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/data/barcode_runs.csv)<br>
List of Illumina barcode samples to be analyzed by the snakemake workflow.<br>
[data/processed_ccs.csv.gz](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/data/processed_ccs.csv.gz)<br>
Processed PacBio CCSs, generated from our [PacBio](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/tree/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/PacBio) routine. Ensure the library is consistent with those used for the assay.<br>
[data/wildtype_sequence.fasta](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/data/wildtype_sequence.fasta)<br>
SARS-CoV-2 nucleocapsid wildtype sequence (Wuhan).<br>

### Sequencing Data

The workflow operates on Illumina barcode sequencing data in fastq.gz format and these files are kept compressed throughout. File location and name should match the listings given in [data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/data/barcode_runs.csv). These files are too large to be contained in GitHub, and so are found, respectively, at:

**p21269-s002_211116_2_S96_L001_R2_001.fastq.gz**<br>
Search NCBI BioProject PRJNA1216977 (BioSample SAMN47131533) or click [here](https://www.ncbi.nlm.nih.gov/bioproject/1216977).

**p21269-s004_211116_4_S98_L001_R2_001.fastq.gz**<br>
Search NCBI BioProject PRJNA1216977 (BioSample SAMN47131889) or click [here](https://www.ncbi.nlm.nih.gov/bioproject/1216977).

## Workflow

Use the `snakemake` environment:

`conda activate snakemake`

Run `snakemake` using specified number of cores:

`snakemake -j 6`

## Key Output
