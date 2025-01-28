# SARS-CoV-2 Nucleocapsid Structure DMS

## Paper
Paper can be found [here]need to provide link when available.

## Analysis Workflow
### Build Computing Environment
`conda` is required, and can be obtained _via_ the minimal installer, `miniconda`, [here](https://docs.anaconda.com/free/miniconda/).

Installation of `mamba` (another package manager) is then required since snakemake depends on it. `mamba` sometimes will not install with default conda settings. A workaround is to change conda settings to 4.12.

`conda install conda=4.12`

`mamba` installation:

`conda install -n base -c conda-forge mamba`

Installation of `snakemake` itself is then required. Ensure that conda is in its base environment:

`conda activate base`

`snakemake` installation:

`mamba create -c conda-forge -c bioconda -n snakemake snakemake`

Now navigate to the new `snakemake` environment:

`conda activate snakemake`

Now that the environment has been set up, various modules need to be installed. Some dependency conflicts can arise. The `dms_variants` and `alignparse` packages were specifically designed for `python 3.7`, but the snakemake package typically installs with a newer version of python, and conflicts can arise if we start to mix different python versions within the one environment. I’ve found that the best way round this issue is to install everything necessary using the default python version which comes with snakemake \(in our case, `python 3.8`\). Another complication is that most of the required packages install with `pip`, not `conda`, but this does not seem to present any issues. The required packages are:

`pip install jupyter`

`pip install dms_variants`

`pip install alignparse`

And, for `seaborn`, it is necessary to install to a specified folder \(i.e. the one `snakemake` will point to\):

`pip install --target=/usr/lib/python3.8 seaborn`

When running `snakemake`, it will most likely throw up an error in the package `dna_features_viewer`, which installs along with `alignparse` \(and is called in some of the Jupyter notebooks, so is a required module!\). The error will most likely present as something like this:

```
File /usr/lib/python3.8/dna_features_viewer/biotools.py:38
     34     return complement(sequence)[::-1]
     36 print(sys.version)
     37 aa_short_to_long_form_dict = {
---> 38     _aa1: _aa3[0] + _aa3[1:].lower() for (_aa1, _aa3) in zip(aa1 + "*", aa3 + ["*"])
     39 }
     42 def translate(dna_sequence, long_form=False):
     43     """Translate the DNA sequence into an amino-acids sequence MLKYQT...
     44 
     45     If long_form is true, a list of 3-letter amino acid representations
     46     is returned instead (['Ala', 'Ser', ...]).
     47     """

TypeError: can only concatenate tuple (not "str") to tuple
TypeError: can only concatenate tuple (not "str") to tuple
```

This error appears in a function which is not necessary for our analysis. Therefore, we can manipulate `biotools.py` (in our case, contained in folder, `/usr/lib/python3.8/dna_features_viewer/` but it could be in a different location depending on where these packages installed). We then comment out the reverse_complement function:

```
# def reverse_complement(sequence):
#    """Return the reverse-complement of the DNA sequence.
#
#    For instance ``complement("ATGCCG")`` returns ``"GCCGTA"``.
#
#    Uses BioPython for speed.
#    """
#    return complement(sequence)[::-1]
#
#aa_short_to_long_form_dict = {
#    _aa1: _aa3[0] + _aa3[1:].lower() for (_aa1, _aa3) in zip(aa1 + "*", aa3 + ["*"])
#}
```
## Twist Library

For the Twist library used, click [here](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/tree/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Twist_Library_QC).

## PacBio Library Sequencing

PacBio library sequencing is required to connect variants to their respective barcodes. This process is described in [PacBio_Library_Sequencing](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/tree/main/PacBio_Library_Sequencing) within the [DMS_IgG1Fc](https://github.com/Ortlund-Laboratory/DMS_IgG1Fc/tree/main) repository.

## Workflow

Input files and scripts for the calculation and presentation of binding scores are provided in [Workflow]need to update this.

## Deposited Data

Relevant input/output data are deposited in the [Deposited_Data]need to add subfolder. Large sequencing files which cannot be maintained on GitHub are provided in external databases, and links are provided. Between these and the example workflows, users should be able to recreate our results. Should issues arise, please contact adkeith@emory.edu or eortlun@emory.edu. 
