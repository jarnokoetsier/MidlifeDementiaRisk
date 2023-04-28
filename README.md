# DNA Methylation-based Prediction of Midlife Dementia Risk
Jarno Koetsier<sup>1</sup>, Rachel Cavill<sup>2</sup>, and Ehsan Pishva<sup>3</sup>

---
![licence](https://badgen.net/badge/Licence/MIT/purple)
![status](https://badgen.net/badge/Status/Complete/green)

This repository contains the scripts that were used for the thesis project *DNA Methylation-based Prediction of Midlife Dementia risk* as part of the Master Systems Biology at Maastricht University.

## Content
1. [Background](#Background)
2. [Research Questions](#Research-Questions)
3. [Workflow](#Workflow)
4. [Contact](#Contact)

## Background
For the early identification of people with increased dementia risk, [LIBRA](https://onlinelibrary.wiley.com/doi/full/10.1002/gps.4245) and [CAIDE](https://www.sciencedirect.com/science/article/pii/S1474442206705373) midlife dementia risk scores have been developed. As DNA methylation may act as the molecular link between lifestyle/environment and the biological processes governing health and disease, DNA methylation data might be utilized instead to predict a midlife dementia risk. However, the large dimensionality of DNA methylation data makes parameter optimization and model training often computationally infeasible without prior feature selection.

![Methylation](/Images/Methylation.PNG?raw=true "Methylation")

## Research Questions
1.	Can a robust and computationally feasible **feature selection** method be established to reduce the dimensionality of the data?
2.	Can reliable DNA methylation-based models for the **prediction** of a person’s LIBRA, CAIDE1, and CAIDE2 scores and risk factors be constructed in a general population cohort?
3.	Does the extension of the dementia risk score models with **polygenetic risk scores (PGSs)** of dementia risk factors, subtypes, and comorbidities improve predictive power? 
4.	Can the LIBRA, CAIDE1, CAIDE2, and risk factor models be used for estimation of a person’s **dementia risk** in dementia-associated cohorts?
5.	What **biological processes** are captured by the most important features of the risk factor models? 

## Workflow
An overview of the applied methodological workflow is shown in figure below and encompasses of the following steps:

### I. Pre-processing
1. Pre-processing of **genomics** data ([GenotypeProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/GenotypeProcessing))
2. Pre-processing of **DNA methylation** data ([MethylationProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/MethylationProcessing))
3. Pre-processing of **phenotype** data ([PhenotypeProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/PhenotypeProcessing)), 

### II. Evaluation of Feature Selection Methods
Evaluation of eight different feature selection methods ([FeatureSelection](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/FeatureSelection)):
1. **Variance**-based feature selection
   * &beta;-values
   * M-values
   * &beta;-values corrected for cell type composition
   * M-values corrected for cell type composition
2. **S-score**-based feature selection
3. **PCA**-based feature selection
4. **Kennard-Stone-like** feature selection
5. **Correlation**-based feature selection

### III. Prediction of Dementia Risk
Prediction of the **dementia risk** scores, categories, and factors ([MachineLearning](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/MachineLearning)). 

### IV. Model Validation & Interpretation
1. **Biological interpretation** of the best-performing risk factor models ([ModelInterpretation](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/ModelInterpretation))
2. **Validation** of the established models in independent dementia-associated cohorts ([ModelValidation](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/ModelValidation)). 


![Workflow](/Images/Workflow.png?raw=true "Workflow")

## Contact
Feel free to contact me via email: j.koetsier@student.maastrichtuniversity.nl

