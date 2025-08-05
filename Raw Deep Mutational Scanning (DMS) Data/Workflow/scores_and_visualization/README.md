# File Prep

[variant_counts.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/results/counts/variant_counts.csv) from results/counts needs to be manipulated.

By referring to [barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/data/barcode_runs.csv), we can see that three experiments were run: a reference sample (exp1-none-0-reference), an enhanced stability sample (exp2-mycpos-2000-escape) and a diminished stability sample (exp3-mycneg-2000-escape). 

**NOTE**: Though the above high/low-Myc samples have escape in their titles, this was due to naming conventions in some of the Jupyter scripts, which requires samples to be classified as either 'reference' or 'escape'.

For compatibility with our R scripts, we require the information in [variant_counts.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/results/counts/variant_counts.csv) to be split into three-column (barcode,mutation,count) files. This is very simple and you can write a script to do this if you like.

First, remove all variants where the number of amino acid mutations is 0 (i.e. is still wildtype) or greater than 1 (multiple point mutations).

```
awk -F',' '!($10!=1)' variant_counts.csv > tmp.csv && mv tmp.csv variant_counts.csv
```
Next step is to separate into reference, mycpos and mycneg files.
```
sed -n '/exp1-none-0-reference/p' variant_counts.csv > ref_variant_counts.txt
sed -n '/exp2-mycpos-2000-escape/p' variant_counts.csv > mycpos_variant_counts.txt
sed -n '/exp3-mycneg-2000-escape/p' variant_counts.csv > mycneg_variant_counts.txt
```
Then reformat each file so that they just contain tab-separated columns for barcode, mutation and count.
```
awk '{print $4, $8, $5}' FS="," OFS="\t" ref_variant_counts.txt > tmp.txt && mv tmp.txt ref_variant_counts.txt
awk '{print $4, $8, $5}' FS="," OFS="\t" mycpos_variant_counts.txt > tmp.txt && mv tmp.txt mycpos_variant_counts.txt
awk '{print $4, $8, $5}' FS="," OFS="\t" mycneg_variant_counts.txt > tmp.txt && mv tmp.txt mycneg_variant_counts.txt
```
Finally, add tab-separated headers (barcode  mutation  count) for each of these files.

Now that we have the input files formatted correctly, we are ready to calculate enhanced/diminished stability scores and to visualize the data. We run these steps in separate mycpos and mycneg subdirectories.

## Z-normalized MycPos minus MycNeg Data

For Figure 2 of our paper, we display our dataset encompassing both stabilization and destabilization scores. To generate this dataset, we z-normalize the MycPos and MycNeg datasets, subtract MycNeg from MycPos, then z-normalize again. This procedure is summarized in [z_norm_mycpos_minus_mycneg.xlsx](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/z_norm_pos_minus_neg/z_norm_mycpos_minus_mycneg.xlsx). See the [z_norm_pos_minus_neg](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/z_norm_pos_minus_neg) subfolder.
