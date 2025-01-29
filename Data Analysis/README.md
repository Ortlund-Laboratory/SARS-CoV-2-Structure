# Data Analysis / Method Normalization & Comparison

## Paper

See Methods in our (INSERT LINK TO PAPER HERE) for details.

## Results

Our results are split into two excel files, for the DD and RBD respectively.

In each file, the first tab **Experiment MycPos - MycNeg** gives the raw score per variant according to MycPos selection in the first two columns. Moving rightwards across the sheet, this data is z-normalized and provided in a separate column. Next, the raw score per variant according to MycNeg selection is given. This too is z-normalized and the new data provided in a separate column. The z-normalized MycNeg data is then subtracted from the z-normalized MycPos data. Finally, this new dataset is itself z-normalized (i.e. a second round of z-normalization is performed on the combined data).

The second tab **ProteinMPNN** gives the raw score per variant in the first two columns. Moving rightwards across the sheet, this data is z-normalized and provided in a separate column. Signs are then reversed so that a higher score indicates a more stable structure and a lower score a less stable structure, to make the data consistent with DMS.

The third tab **Rosetta** gives the raw score per variant in the first two columns. Moving rightwards across the sheet, this data is z-normalized and provided in a separate column. Signs are then reversed so that a higher score indicates a more stable structure and a lower score a less stable structure, to make the data consistent with DMS.

The fourth tab **All** combines the data from the first three tabs. Care is taken that any rows where there is missing data (which occasionally occurs when a DMS variant is not detected) are deleted. Average scores according to all three methods are calculated per variant. Per-site scores for DMS, ProteinMPNN, Rosetta and an average of all three are also calculated.

The fifth tab **Difference DMS Rosetta** calculates the difference in predicted Rosetta score and generated DMS score per variant (DMS minus Rosetta). These differences are then sorted from most negative (Rosetta predicts mutation to be stabilizing but DMS finds it to be destabilizing) to most positive (*vice versa*) values. 
