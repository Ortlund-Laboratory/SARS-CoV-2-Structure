# Data Analysis / Method Normalization & Comparison

## Paper

See Methods in our (INSERT LINK TO PAPER HERE) for details.

## Results

Our results are split into two excel files, for the DD and RBD respectively.

In each file, the first tab **Experiment MycPos - MycNeg** gives the raw score per variant according to MycPos selection in the first two columns. Moving rightwards across the sheet, this data is z-normalized and provided in a separate column. Next, the raw score per variant according to MycNeg selection is given. This too is z-normalized and the new data provided in a separate column. The z-normalized MycNeg data is then subtracted from the z-normalized MycPos data. Finally, this new dataset is itself z-normalized (i.e. a second round of z-normalization is performed on the combined data).

The second tab **ProteinMPNN** gives the raw score per variant in the first two columns. 
