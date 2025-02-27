## Input Files Required

[gfp_myc_analysis_BarcodeMapping.R](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/scores_and_visualization/mycpos/RBD_only/gfp_myc_analysis_BarcodeMapping.R)<br>
R script to calculate scores from reference and mycpos counts, generate heatmaps and produce files to map scores onto structures.<br>
[ref_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/scores_and_visualization/mycpos/RBD_only/ref_variant_counts.txt)<br>
Barcode counts for the reference sample.<br>
[mycpos_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/scores_and_visualization/mycpos/RBD_only/mycpos_variant_counts.txt)<br>
Barcode counts for the mycpos sample.<br>
[N_Wuhan.fasta](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/scores_and_visualization/mycpos/RBD_only/N_Wuhan.fasta)<br>
Amino acid sequence for N Wuhan. This is required to complete the heatmap.

## Workflow

```
rstudio gfp_myc_analysis_BarcodeMapping.R
```

## Key Output

[GFP+_Myc+_scores.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/scores_and_visualization/mycpos/output/GFP+_Myc+_scores.csv)<br>
Log of each mutation and its associated score.<br>
[GFP+_Myc+_per_site_score.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/scores_and_visualization/mycpos/output/GFP+_Myc+_per_site_score.csv)<br>
Log of each site and its associated averaged score.<br>
[GFP+_Myc+_heatmap01.png](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/scores_and_visualization/mycpos/output/GFP+_Myc+_heatmap01.png)<br>
Heatmap of scores for all single-point variants from the MycPos sample (sites 2-208).<br>
[GFP+_Myc+_heatmap02.png](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Initial_GFP+_Myc+_Selection/scores_and_visualization/mycpos/output/GFP+_Myc+_heatmap02.png)<br>
Heatmap of scores for all single-point variants from the MycPos sample (sites 211-419).<br>

