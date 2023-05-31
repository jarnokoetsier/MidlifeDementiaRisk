library(tidyverse)
library(caret)


# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("~/allModels.RData")
load("PPMI/metaData_ppmi.RData")
load("PPMI/X_PPMI_imp.RData")

# Predict factors
factors <- names(allModels)
predictedScore_factors <- matrix(NA, nrow = nrow(X_PPMI_imp), ncol = length(factors))
colnames(predictedScore_factors) <- factors
rownames(predictedScore_factors) <- rownames(X_PPMI_imp)
for (f in 1:length(factors)){
  f1 <- factors[f]
  model <- allModels[[f1]]
  
  # Get features for model fitting
  features <- colnames(model$trainingData)[-10001]
  
  # Make predictions
  predictedScore_factors[,f] <- predict(model, X_PPMI_imp[,features], type = "prob")$No
  
}
predictedScore_factors <- as.data.frame(predictedScore_factors)


# Predict age
X_PPMI_imp_B <- (2^X_PPMI_imp)/(1+ 2^X_PPMI_imp)
predictedAge <- agep(t(X_PPMI_imp_B), 
                     method='all')

predictedScore_factors$Age <- predictedAge$skinblood.skinblood.age

# Save predicted scores
save(predictedScore_factors, file = "PPMI/predictedScore_factors_PPMI.RData")



                     