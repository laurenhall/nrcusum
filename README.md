# nrcusum
Supporting materials for *A modified CUSUM test to control post-outbreak false alarms.*

## Files included
**nrcusum_functions.R** : R code containing two functions for implementing the corrected CUSUM. Requires the `pbapply` package.
**nrcusum_correction.R** : R code for generation of simulated data and implementation of CUSUM correction. Includes the generation of figures for main text. Requires the functions from `nrcusum_functions.R` and the `pROC` package for AUC analysis.
**germany_salmonella_example.R** : R code for main text data demonsration. Requires the functions from `nrcusum_functions.R` and the `surveillance` package, which includes the data. 
**simdata.rda** : The data sets simulated by `nrcsum_correction.R`. 
