# WM-connectome-devlopment
This repository provides data and code for reproducing a range of analyses involving our study.
If you use these data and code, please cite our paper:
Guozheng Feng, Rui Chen, Rui Zhao, et al. Longitudinal development of the human white matter structural connectome and its association with brain transcriptomic and cellular architecture. Communications Biology.

https://doi.org/10.5281/zenodo.10212534

### Dependencies:
The lmer function in R4.1.2 software (https://www.r-project.org/) is used to perform mixed effect model. 
Matlab scripts to run the preprocessing of AHBA dataset can be found at https://github.com/BMHLab/AHBAprocessing. 

### Data and code files
1. 'preprocess.txt' - A script is used to preprocess MRI.
2. 'MLM.R' - Estimating both linear and quadratic models by a mixed effect model.
3. 'slopeGenePLS.m' - A script is used to perform PLS correlation.
4. 'metascape_result.xlsx' - Gene functional enrichment results for the GO biological process pathway search with Metascape
