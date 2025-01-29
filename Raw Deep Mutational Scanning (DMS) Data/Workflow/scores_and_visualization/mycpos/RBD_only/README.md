## Input Files Required

[RBD_only_mycpos.R](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycpos/RBD_only/RBD_only_mycpos.R)<br>
R script to calculate scores from reference and mycpos counts, generate heatmaps and produce files to map scores onto structures.<br>
[ref_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycpos/RBD_only/ref_variant_counts.txt)<br>
Barcode counts for the reference sample.<br>
[mycpos_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycpos/RBD_only/mycpos_variant_counts.txt)<br>
Barcode counts for the mycpos sample.<br>
[N_Wuhan.fasta](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycpos/RBD_only/N_Wuhan.fasta)<br>
Amino acid sequence for N Wuhan. This is required to complete the heatmap.

## Workflow

```
rstudio RBD_only_mycpos.R
```

## Key Output

[Mycpos_scores.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycpos/RBD_only/output/Mycpos_scores.csv)<br>
Log of each mutation and its associated score.<br>
[Mycpos_per_site_score.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycpos/RBD_only/output/Mycpos_per_site_score.csv)<br>
Log of each site and its associated averaged score.<br>
[Mycpos_Fraction_heatmap01.png](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycpos/RBD_only/output/Mycpos_Fraction_heatmap01.png)<br>
Heatmap of scores for all single-point variants from the MycPos sample.<br>
