# WM-connectome-devlopment
This repository provides data and code for reproducing a range of analyses involving our study.
If you use these data and code, please cite our paper:
Guozheng Feng, Rui Chen, Rui Zhao, et al. Longitudinal development of the human white matter structural connectome and its association with brain transcriptomic and cellular architecture. Communications Biology.

### Dependencies:
The lmer function in R4.1.2 software (https://www.r-project.org/) is used to perform mixed effect model. 
Matlab scripts to run the preprocessing of AHBA dataset can be found at https://github.com/BMHLab/AHBAprocessing. 

### Data and code files
1. 'MLM.R' - Estimating both linear and quadratic models by a mixed effect model.
2. 'fnfa246_6-13slope.mat' - developmental slopes of regional network properties (nodal efficiency, nodal local efficiency and nodal degree centrality)
3. 'ge_slope_surrogates.mat' - 10,000 surrogate maps of developmental slope of nodal efficiency by a spatially constrained generation model
4. 'AnalysisSlope.m' - A script is used to analyze the significance in different properties and the relationship with AHBA data
5. 'metascape_result.xlsx' - Gene functional enrichment results for the GO biological process pathway search with Metascape
