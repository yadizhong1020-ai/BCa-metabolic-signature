Overview

This repository contains the R and Python scripts used for the analyses presented in our study investigating longitudinal metabolic alterations and metabolomic classification in non-muscle-invasive bladder cancer (NMIBC).

The analytical workflow includes:

- UK Biobank longitudinal metabolomic trajectory analysis
- 1:5 case–control matching
- Covariate adjustment and residual calculation
- Time-resolved metabolite trajectory analysis
- Metabolite trajectory clustering
- DE-SWAN analysis
- Consensus feature selection (LASSO + SVM-RFE + Boruta)
- Random Forest-based metabolomic classification model
- Bootstrap evaluation of model performance

Software Requirements
R
R >= 4.3

Packages:
- data.table
- dplyr
- caret
- glmnet
- Boruta
- randomForest
- pROC
- ggplot2
- pheatmap
Python
Python >= 3.10

Packages:
- numpy
- pandas
- scipy
- scikit-learn
- matplotlib
- seaborn
- lifelines
