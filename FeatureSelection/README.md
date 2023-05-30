# Feature selection
* `CorVsNon_Features.R`: Compare feature selection by ElasticNet with and without correlation-based feature selection.
* `EvaluateFeatureSelection.R`: Evaluate different feature selection methods in terms of common features, overlapping information, and explained variance by cell type composition.
* `EvaluateFeatureSelection_CAIDE1.R`: Evaluate different feature selection methods in terms of CAIDE1 predictive performance. For this, `MachineLearning/Score/PredictCAIDE1_Score.R` script has to be run first.
* `EvaluateGreedy.R`: Evaluate Greedy feature selection.
* `FeatureSelection_Correlation.R`: Perform correlation-based feature selection.
* `FeatureSelection_Unsupervised.R`: Perform all unsupervised feature selection methods: variance, S-score, KS-like, PCA.
* `GOenrichment.R`: Perform GO overrepresentation analysis on the selected features.
* `Greedy.R`: Perform Greedy feature selection.
