See [z_norm_mycpos_minus_mycneg.xlsx](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/z_norm_pos_minus_neg/z_norm_mycpos_minus_mycneg.xlsx) for the z-normalization process. After this, We use the following input files:

## Input Files Required

[MutationMapping_integrated_scores_DMS_enrich_minus_escape.R](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/z_norm_pos_minus_neg/MutationMapping_integrated_scores_DMS_enrich_minus_escape.R)<br>
R script to plot MycPos-MycNeg scores.<br>
[MycPos_Minus_MycNeg_Scores.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/z_norm_pos_minus_neg/MycPos_Minus_MycNeg_Scores.csv)<br>
Z-normalized mycpos-myc neg scores.<br>
[N_Wuhan.fasta](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/z_norm_pos_minus_neg/N_Wuhan.fasta)<br>
Amino acid sequence for N Wuhan. This is required to complete the heatmap.

## Workflow

```
rstudio MutationMapping_integrated_scores_DMS_enrich_minus_escape.R
```

## Key Output

[Integrated_Nucleocapsid_DMS_mycpos_minus_mycneg_per_site_score.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/z_norm_pos_minus_neg/Integrated_Nucleocapsid_DMS_mycpos_minus_mycneg_per_site_score.csv)<br>
Log of each site and its associated averaged score.<br>
[Integrated_Nucleocapsid_DMS_mycpos_minus_mycneg_heatmap01.png](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/z_norm_pos_minus_neg/Integrated_Nucleocapsid_DMS_mycpos_minus_mycneg_heatmap01.png)<br>
Heatmap of scores for all single-point variants from the MycPos-MycNeg sample (part 1).<br>
[Integrated_Nucleocapsid_DMS_mycpos_minus_mycneg_heatmap02.png](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/z_norm_pos_minus_neg/Integrated_Nucleocapsid_DMS_mycpos_minus_mycneg_heatmap02.png)<br>
Heatmap of scores for all single-point variants from the MycPos-MycNeg sample (part 2).<br>
