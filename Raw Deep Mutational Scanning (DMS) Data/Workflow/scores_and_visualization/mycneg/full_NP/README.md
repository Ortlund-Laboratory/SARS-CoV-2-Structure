## Input Files Required

[Mar2022_mycneg_myBarcodeMapping.R](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/Mar2022_mycneg_myBarcodeMapping.R)<br>
R script to calculate scores from reference and mycneg counts, generate heatmaps and produce files to map scores onto structures.<br>
[ref_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/ref_variant_counts.txt)<br>
Barcode counts for the reference sample.<br>
[mycneg_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/mycneg_variant_counts.txt)<br>
Barcode counts for the mycneg sample.<br>
[N_Wuhan.fasta](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/N_Wuhan.fasta)<br>
Amino acid sequence for N Wuhan. This is required to complete the heatmap.

## Workflow

```
rstudio Mar2022_mycneg_myBarcodeMapping.R
```

## Key Output

[MycNeg_scores.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/output/MycNeg_scores.csv)<br>
Log of each mutation and its associated score.<br>
[MycNeg_per_site_score.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/output/MycNeg_per_site_score.csv)<br>
Log of each site and its associated averaged score.<br>
[MycNeg_Fraction_heatmap01.png](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/output/MycNeg_Fraction_heatmap01.png)<br>
Heatmap of scores for all single-point variants from the MycNeg sample (part 1).<br>
[MycNeg_Fraction_heatmap02.png](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/output/MycNeg_Fraction_heatmap02.png)<br>
Heatmap of scores for all single-point variants from the MycNeg sample (part 2).<br>
