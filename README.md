![licence](https://badgen.net/badge/Licence/MIT/purple)
![status](https://badgen.net/badge/Status/Complete/green)

<h1 align="center">
DNA Methylation-based Prediction
   <br>
of Midlife Dementia Risk
</h1>

<p align="center">
<a href="https://github.com/jarnokoetsier/MidlifeDementiaRisk/blob/main/README.md#Background">Background</a>
     ·
<a href="https://github.com/jarnokoetsier/MidlifeDementiaRisk/blob/main/README.md#Research-Questions">Research Questions</a>
     ·
<a href="https://github.com/jarnokoetsier/MidlifeDementiaRisk/blob/main/README.md#Workflow">Workflow</a>
     ·
<a href="https://github.com/jarnokoetsier/MidlifeDementiaRisk/blob/main/README.md#Contributors">Contributors</a>
</p>

<p align="center">
This repository contains the scripts that were used for the thesis project <i>DNA Methylation-based Prediction of Midlife Dementia risk</i> as part of the Master Systems Biology at Maastricht University.
</p>

 <br>

## Background
For the early identification of people with increased dementia risk, [LIBRA](https://onlinelibrary.wiley.com/doi/full/10.1002/gps.4245) and [CAIDE](https://www.sciencedirect.com/science/article/pii/S1474442206705373) midlife dementia risk scores have been developed. As DNA methylation may act as the molecular link between lifestyle/environment and the biological processes governing health and disease, DNA methylation data might be utilized instead to predict a midlife dementia risk. However, the large dimensionality of DNA methylation data makes parameter optimization and model training often computationally infeasible without prior feature selection.

![Methylation](/Images/Methylation.PNG?raw=true "Methylation")

## Research Questions
1.	Can a robust and computationally feasible **feature selection** method be established to reduce the dimensionality of the data?
2.	Can reliable DNA methylation-based models for the **prediction** of a person’s LIBRA, CAIDE1, and CAIDE2 scores and risk factors be constructed in a general population cohort?
3.	Does the extension of the dementia risk score models with **polygenetic risk scores (PGSs)** of dementia risk factors, subtypes, and comorbidities improve predictive power? 
4.	Can the LIBRA, CAIDE1, CAIDE2, and risk factor models be used to predict **dementia** and **cognitive impairment** status in dementia-associated cohorts??
5.	What **biological processes** are captured by the most important features of the risk factor models? 

## Workflow
An overview of the applied methodological workflow is shown in figure below and encompasses of the following steps:

### I. Pre-processing
1. Pre-processing of **phenotype** data ([PhenotypeProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/1.%20PhenotypeProcessing)) 
2. Pre-processing of **genomics** data ([GenotypeProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/2.%20GenotypeProcessing))
3. Pre-processing of **DNA methylation** data ([MethylationProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/3.%20MethylationProcessing))


### II. Evaluation of Feature Selection Methods
Evaluation of eight different feature selection methods ([FeatureSelection](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/4.%20FeatureSelection)):
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
Prediction of the **dementia risk** scores, categories, and factors ([MachineLearning](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/5.%20MachineLearning)). 

### IV. Model Validation & Interpretation
1. **Biological interpretation** of the best-performing risk factor models ([ModelInterpretation](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/6.%20ModelInterpretation))
2. **Validation** of the established models in independent dementia-associated cohorts ([ModelValidation](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/7.%20ModelValidation)). 


![Workflow](/Images/Workflow.png?raw=true "Workflow")

## Contributors
Jarno Koetsier<sup>1*</sup>, Rachel Cavill<sup>2</sup>, and Ehsan Pishva<sup>3</sup>

*<sup>1</sup> Faculty of Science and Engineering (FSE), Maastricht University*

*<sup>2</sup> Department of Advanced Computing Sciences (DACS), Faculty of Science and Engineering (FSE), Maastricht University*

*<sup>3</sup> Department of Psychiatry and Neuropsychology, School for Mental Health and Neuroscience (MHeNS), Faculty of Health, Medicine and Life Sciences (FHML), Maastricht University*

---

<sup>*</sup>Feel free to contact me via email: j.koetsier@student.maastrichtuniversity.nl

