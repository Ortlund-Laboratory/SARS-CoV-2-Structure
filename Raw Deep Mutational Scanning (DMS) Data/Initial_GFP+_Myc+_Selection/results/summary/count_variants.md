# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Initialize codon variant table
Initialize the [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) using the wildtype gene sequence and the CSV file with the table of variants:


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons'],
                                     feature_parse_specs=config['feature_parse_specs'])
geneseq = targets.get_target(config['primary_target']).get_feature('gene').seq
print(f"Read gene of {len(geneseq)} nt for {config['primary_target']} from {config['amplicons']}")
      
print('Initializing CodonVariantTable from gene sequence and ' +
      config['codon_variant_table'])
      
variants = dms_variants.codonvarianttable.CodonVariantTable(
                geneseq=geneseq,
                barcode_variant_file=config['codon_variant_table'],
                substitutions_are_codon=True,
                substitutions_col='codon_substitutions',
                primary_target=config['primary_target'])
```

    Read gene of 1257 nt for nucleocapsid from data/PacBio_amplicons.gb
    Initializing CodonVariantTable from gene sequence and results/prior_DMS_data/codon_variant_table.csv


## Setup to parse barcodes
Read data frame with list of all barcode runs.
Note how multiple R1 files are delimited by `; ` and are split out separately:


```python
print(f"Reading list of barcode runs from {config['rt_barcode_runs']}")

rt_barcode_runs = (pd.read_csv(config['rt_barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )
      
display(HTML(rt_barcode_runs.to_html(index=False)))
```

    Reading list of barcode runs from data/rt_barcode_runs.csv



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>date</th>
      <th>experiment</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>concentration_units</th>
      <th>group</th>
      <th>selection</th>
      <th>frac_escape</th>
      <th>cells_sorted</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
      <td>250224</td>
      <td>exp1</td>
      <td>none</td>
      <td>0</td>
      <td>ng_per_mL</td>
      <td>none</td>
      <td>reference</td>
      <td>NaN</td>
      <td>5000000</td>
      <td>[/home/adkeith@Eu.Emory.Edu/DMS_Workflow/nucleocapsid/p21269_gfp_myc_analysis/rt_p21269-s002_211116_2_S96_L001_R2_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>N-Lib05</td>
      <td>exp2-mycpos-2000-escape</td>
      <td>250224</td>
      <td>exp2</td>
      <td>mycpos</td>
      <td>2000</td>
      <td>ng_per_mL</td>
      <td>DK</td>
      <td>escape</td>
      <td>0.15</td>
      <td>1000000</td>
      <td>[/home/adkeith@Eu.Emory.Edu/DMS_Workflow/nucleocapsid/p21269_gfp_myc_analysis/rt_p21269-s004_211116_4_S98_L001_R2_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(rt_barcode_runs) == len(rt_barcode_runs.groupby(['library', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(rt_barcode_runs['library']) - set(variants.libraries)
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGAGAGGGGCGGNATCCGCGGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td>TTAACGTGGCTTCTTCTGCCACAGCATGATGAGAATAATAAGGGAAATGATAG</td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>10</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>10</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **parser_params)
           for lib in variants.libraries}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>N-Lib05</td>
      <td>156624</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 12 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.sample)
                  for run in rt_barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>CCCCGCCTTGGACCG</td>
      <td>3280</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
    <tr>
      <td>CGGTGCGCTCCGGGC</td>
      <td>3078</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
    <tr>
      <td>CGGCGTGGCACGCTC</td>
      <td>2830</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
    <tr>
      <td>TGACAGGATGAGCCG</td>
      <td>2821</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
    <tr>
      <td>CCGCCCCATTTGCAG</td>
      <td>2738</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>45453264</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>10943262</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>7648701</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>5358820</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">N-Lib05</th>
      <th>exp1-none-0-reference</th>
      <td>0</td>
      <td>7648701</td>
      <td>10943262</td>
      <td>5358820</td>
      <td>45453264</td>
    </tr>
    <tr>
      <th>exp2-mycpos-2000-escape</th>
      <td>0</td>
      <td>1757477</td>
      <td>2595861</td>
      <td>1449388</td>
      <td>11157227</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library') +
    facet_grid('sample ~ library') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_36_0.png)
    


## Add barcode counts to variant table
Now we use the [CodonVariantTable.add_sample_counts_df](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable.add_sample_counts_df) method to add the barcode counts to the variant table:


```python
variants.add_sample_counts_df(counts)
```

The variant table now has a `variant_count_df` attribute that gives a data frame of all the variant counts.
Here are the first few lines:


```python
display(HTML(variants.variant_count_df.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>sample</th>
      <th>barcode</th>
      <th>count</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>nucleocapsid</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
      <td>CCCCGCCTTGGACCG</td>
      <td>3280</td>
      <td>33</td>
      <td>CCA20CTG TTT110GTT GAG378TGA</td>
      <td>P20L F110V E378*</td>
      <td>3</td>
      <td>3</td>
    </tr>
    <tr>
      <td>nucleocapsid</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
      <td>CGGTGCGCTCCGGGC</td>
      <td>3078</td>
      <td>34</td>
      <td>GAG31TTC</td>
      <td>E31F</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>nucleocapsid</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
      <td>CGGCGTGGCACGCTC</td>
      <td>2830</td>
      <td>23</td>
      <td>GAG323ATG</td>
      <td>E323M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>nucleocapsid</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
      <td>TGACAGGATGAGCCG</td>
      <td>2821</td>
      <td>24</td>
      <td>GAT288GTG</td>
      <td>D288V</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>nucleocapsid</td>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
      <td>CCGCCCCATTTGCAG</td>
      <td>2738</td>
      <td>21</td>
      <td>TTT315ACC</td>
      <td>F315T</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


Confirm that we have counts for all of our library / sample combinations:


```python
pd.testing.assert_frame_equal(
    rt_barcode_runs[['library', 'sample']]
                .sort_values(['library', 'sample'])
                .reset_index(drop=True),
    variants.variant_count_df
        [['library', 'sample']]
        .drop_duplicates()
        .sort_values(['library', 'sample'])
        .reset_index(drop=True),
    check_dtype=False,
    check_categorical=False,
    check_like=True,
    )
```


```python
rt_barcode_runs[['library', 'sample']].sort_values(['library', 'sample'])
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>N-Lib05</td>
      <td>exp1-none-0-reference</td>
    </tr>
    <tr>
      <th>1</th>
      <td>N-Lib05</td>
      <td>exp2-mycpos-2000-escape</td>
    </tr>
  </tbody>
</table>
</div>




```python
((variants.variant_count_df
        [['library', 'sample']]
        .drop_duplicates()
        .sort_values(['library', 'sample'])
        .reset_index(drop=True)) == rt_barcode_runs[['library', 'sample']].sort_values(['library', 'sample']).reset_index(drop=True)).query('sample == False')
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>
</div>



Write the variant counts data frame to a CSV file.
It can then be used to re-initialize a [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) via its [from_variant_count_df](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable.from_variant_count_df) method:


```python
print(f"Writing variant counts to {config['variant_counts']}")
variants.variant_count_df.to_csv(config['variant_counts'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```


```python

```
