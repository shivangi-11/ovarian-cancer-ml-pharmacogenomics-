# ovarian-cancer-ml-pharmacogenomics-
# Pharmacogenomics-Based Drug Response Prediction in Ovarian Cancer

## Overview
This repository provides the computational workflow used for drug response prediction and biomarker identification in ovarian cancer using pharmacogenomics datasets.

The study integrates CCLE and GDSC datasets, applies machine learning models, and validates biomarkers using TCGA-OV data.

---
 Workflow Summary

The study is divided into three main components:

1. **Data Integration & Consistency Analysis**
   - Integration of CCLE and GDSC datasets
   - Drug sensitivity modeling using AAC
   - Concordance analysis using CI and mCI

2. **Machine Learning Modeling**
   - Preprocessing and batch correction (ComBat)
   - Feature selection using mRMR
   - Regression models (Elastic Net, Ridge, Lasso)

3. **Biomarker Validation**
   - Validation using TCGA-OV dataset
   - Expression analysis and correlation assessment

---

## 📁 Repository Structure
│
├── Scripts/
│ ├── 01_data_concordance_analysis.R
│ ├── 02_ml_preprocessing.R
│ ├── 03_feature_selection.R
│ ├── 04_ml_model_training.R
│ ├── 05_biomarker_identification.R  
│ ├── 06_biomarker_validation.R
├── Results/
│ ├── 01_data_concordance_analysis
│ ├── 02_ml_preprocessing.R
│ ├── 03_feature_selection
│ ├── 04_ml_model_training.R
│ ├── 05_biomarker_identification
│ ├── 06_biomarker_validation
├── README.md

Data Availability

The datasets used in this study are publicly available:

CCLE: https://portals.broadinstitute.org/ccle
GDSC: https://www.cancerrxgene.org/
TCGA-OV: accessed via TCGAbiolinks (For external validation)

** Due to size and licensing restrictions, raw datasets are not included in this repository.
** Also due to size constraints, only processed results for a single drug (Erlotinib) have been included as a representative sample. Additionally, PCA results are also provided as sample outputs.

License

This project is provided for academic and research purposes.
