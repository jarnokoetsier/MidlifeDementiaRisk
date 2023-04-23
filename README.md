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
An overview of the applied methodological workflow is shown in figure below and includes the pre-processing of genomics data ([GenotypeProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/GenotypeProcessing)), DNA methylation data ([MethylationProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/MethylationProcessing)), and phenotype data ([PhenotypeProcessing](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/PhenotypeProcessing)), the evaluation of feature selection methods ([FeatureSelection](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/FeatureSelection)), the prediction of the dementia risk scores and factors ([MachineLearning](https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/MachineLearning)), as well as the biological interpretation ([ModelInterpretation]((https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/ModelInterpretation)) and validation of the established models in independent dementia-associated cohorts ([ModelValidation]((https://github.com/jarnokoetsier/MidlifeDementiaRisk/tree/main/ModelValidation)). The analyses were performed in the data science research environment (DSRI) hosted at Maastricht University.

![Workflow](/Images/Workflow.png?raw=true "Workflow")

## Contact
Feel free to contact me via email: j.koetsier@student.maastrichtuniversity.nl

