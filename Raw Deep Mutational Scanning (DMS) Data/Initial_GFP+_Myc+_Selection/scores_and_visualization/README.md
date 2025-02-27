# File Prep

[variant_counts.csv.gz](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/results/counts/variant_counts.csv.gz) from results/counts needs to be manipulated.

First, gunzip the file, Then, by referring to [barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/data/barcode_runs.csv), we can see that two experiments were run: a reference GFP<sup>+</sup> sample (exp1-none-0-reference) and a GFP<sup>+</sup>Myc<sup>+</sup> sample (exp2-mycpos-2000-escape). 

**NOTE**: Though the above sample has escape in its title, this was due to naming conventions in some of the Jupyter scripts, which requires samples to be classified as either 'reference' or 'escape'.

For compatibility with our R scripts, we require the information in [variant_counts.csv.gz](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/results/counts/variant_counts.csv.gz) to be split into three-column (barcode,mutation,count) files. This is very simple and you can write a script to do this if you like.

First, remove all variants where the number of amino acid mutations is 0 (i.e. is still wildtype) or greater than 1 (multiple point mutations).

```
awk -F',' '!($10!=1)' variant_counts.csv > tmp.csv && mv tmp.csv variant_counts.csv
```
Next step is to separate into reference and sample files
```
sed -n '/exp1-none-0-reference/p' variant_counts.csv > ref_variant_counts.txt
sed -n '/exp2-mycpos-2000-escape/p' variant_counts.csv > mycpos_variant_counts.txt
```
Then reformat each file so that they just contain tab-separated columns for barcode, mutation and count.
```
awk '{print $4, $8, $5}' FS="," OFS="\t" ref_variant_counts.txt > tmp.txt && mv tmp.txt ref_variant_counts.txt
awk '{print $4, $8, $5}' FS="," OFS="\t" mycpos_variant_counts.txt > tmp.txt && mv tmp.txt mycpos_variant_counts.txt
```
Finally, add tab-separated headers (barcode  mutation  count) for each of these files.

Now that we have the input files formatted correctly, we are ready to calculate enhanced expression scores and to visualize the data. We run these steps in a separate mycpos subdirectory.

