# nrcusum
Supporting materials for *A modified CUSUM test to control post-outbreak false alarms.*

## Files included
**nrcusum_functions.R** : R code containing two functions for implementing the corrected CUSUM. Requires the `pbapply` and `pROC` R packages.

**nrcusum_correction.R** : R code for generation of simulated data and implementation of CUSUM correction. Includes the generation of figures for main text. Requires the functions from `nrcusum_functions.R`.

**nrcusum_correction_expanded.R** : An expanded version of `nrcusum_correction.R` which implements each correction separately rather than jointly, making it easier to explore the effects of individual corrections. 

**germany_salmonella_example.R** : R code for main text data demonstration. Requires the functions from `nrcusum_functions.R` and the `surveillance` package, which includes the data. 

**simdata.rda** : The data sets simulated by `nrcsum_correction.R`. 

**cusum_example.rda** : The data set used for the CUSUM example in Figure 1. 
