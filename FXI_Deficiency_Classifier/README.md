# FXI Deficiency Logistic Regression Classifier

This repository contains code and analysis for predicting bleeding phenotype in Factor XI (FXI) deficiency. 

Title: Plasma proteomic signatures determine clot stability and predict bleeding risk in factor XI deficiency

Authors: Stéphanie E. Reitsma, Tracey G. Oellerich, Olivia C. White, Michael L. Hill, Daniel Kissau, John Burthem, Dougald M. Monroe, Paula H. B. Bolton-Maggs, Paul S. de Vries, Suzanne S. Sindi, Karin Leiderman, Gillian N. Gidley, Kirk C. Hansen, Alisa S. Wolberg

## Overview
This repository allows users to build and evaluate logistic regression models to predict FXI deficiency. Key features include:

- **User-defined inputs:**  
  - Fixed features to include in every model  
  - Number of features per model  
  - Output file names  

- **Model evaluation:**  
  - Systematically evaluates all feature combinations using 5-fold cross-validation  
  - Computes performance metrics:
    - Mean AUROC with 95% confidence intervals  
    - Optimal classification thresholds (Youden’s J statistic)  
    - Confusion matrices  
    - Standard rates: TPR, FNR, FPR, TNR  

- **Outputs:**  
  - Excel file containing results  
  - Top 5 models sorted by mean AUROC  
  - ROC curve plots for Manchester and Indeterminate datasets
  
- **Data files included:**  
  - Discovery data: `FXI_Deficient_Bleeders_NonBleeders_data.xlsx`  
  - Hold-out data: `Indeterminant_data.xlsx`

## Installation
Clone the repository and install dependencies:

```bash
git clone https://github.com/karinleiderman/LeidermanResearchGroupCode/tree/main/FXI_Deficiency_Classifier.git
cd FXI_Deficiency_Classifier
pip install -r requirements.txt
```

## Contact
If you have any questions or GitHub issues, feel free to contact us through email (toelleri@unc.edu & karin.leiderman@unc.edu).
