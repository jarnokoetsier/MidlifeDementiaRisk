# DNA Methylation-based Prediction of Midlife Dementia Risk
![licence](https://badgen.net/badge/Licence/MIT/purple)
![status](https://badgen.net/badge/Status/Complete/green)

## Content
1. [Background](#Data)
2. [Research Questions](#Research-aim)
3. [Workflow](#Analysing-the-data)
4. [Contact](#Contact)

## Background

## Research Questions
1.	Can a robust and computationally feasible feature selection method be established to reduce the dimensionality of the data?
2.	Can reliable DNA methylation-based models for the prediction of a person’s LIBRA, CAIDE1, and CAIDE2 scores and risk factors be constructed in a general population cohort?
3.	Does the extension of the dementia risk score models with polygenetic risk scores (PGSs) of dementia risk factors, subtypes, and comorbidities improve predictive power? 
4.	Can the LIBRA, CAIDE1, CAIDE2, and risk factor models be used for estimation of a person’s dementia risk in dementia-associated cohorts?
5.	What biological processes are captured by the most important features of the risk factor models? 

## Workflow
An overview of the applied methodological workflow is shown in figure below and encompasses of the following steps:

### I. Pre-processing
1. Pre-processing of genomics data ([GenotypeProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/GenotypeProcessing))
2. Pre-processing of DNA methylation data ([MethylationProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/MethylationProcessing))
3. Pre-processing of phenotype data ([PhenotypeProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/PhenotypeProcessing)), 

### II. Evaluation of Feature Selection Methods
Evaluation of different feature selection methods ([FeatureSelection](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/FeatureSelection)):
1. Variance-based feature selection
  1. &beta;-values
  2. M-values
  3. &beta;-values + cell type composition correction
  4. M-values + cell type composition correction
2. S-score-based feature selection
3. PCA-based feature selection
4. Kennard-Stone-like feature selection
5. Correlation-based feature selection

### III. Prediction of Dementia Risk
Prediction of the dementia risk scores, categories, and factors ([MachineLearning](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/MachineLearning)). 

### IV. Model Validation & Interpretation
1. Biological interpretation ([ModelInterpretation](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/ModelInterpretation))
2. Validation of the established models in independent dementia-associated cohorts ([ModelValidation](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/ModelValidation)). 


![Workflow](/Images/Workflow.png?raw=true "Workflow")

## Contact
Feel free to contact me via email: j.koetsier@student.maastrichtuniversity.nl

