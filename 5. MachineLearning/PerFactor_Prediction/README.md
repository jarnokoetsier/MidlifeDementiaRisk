# Dementia Risk Factors Predictions
* `/CombineFactors/`: Combine the risk factors for CAIDE1, CAIDE2, and LIBRA prediction
* `AgeAndSmokingPrediction.R`: Evaluate age and smoking models.
* `EvaluatePerFactor.R`: Evaluate all risk factor models.
* `GetBestModels.R`: Get best-performing risk factor models.
* `HillaryModels.csv`: Risk factor models from Hillary & Marioni (2020), doi:10.12688/wellcomeopenres.16458.2.
* `HillaryModels.R`: Construct risk factor models from Hillary & Marioni (2020), 10.12688/wellcomeopenres.16458.2.
* `PerFactorPrediction_EN.R`: Construction of risk factor models by correlation-based feature selection and ElasticNet algorithm.
* `PerFactorPrediction_lit.R`: Construction of risk factor models by literature-based feature selection and ElasticNet algorithm.
* `PerFactorPrediction_RF.R`: Construction of risk factor models by correlation-based feature selection and Random Forest algorithm.
* `ROC_PerFactor.R`: Construction of ROC curves for each risk factor prediction in test and cross-validation.
* `selectedProbes_lit.RData`: Selected probes by literature-based feature selection (genome-wide significance in prior EWASs).
* `SmokingScoreRefData.rda`: Smoking Score model by Elliot et al. (2014), https://doi.org/10.1186/1868-7083-6-4.
