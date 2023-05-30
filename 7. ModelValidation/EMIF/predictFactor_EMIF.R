library(tidyverse)
library(caret)


# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("~/allModels.RData")
load("EMIF/metaData_EMIF.RData")
load("EMIF/X_EMIF_imp.RData")

# Predict factors
factors <- names(allModels)
predictedScore_factors <- matrix(NA, nrow = nrow(X_EMIF_imp), ncol = length(factors))
colnames(predictedScore_factors) <- factors
rownames(predictedScore_factors) <- rownames(X_EMIF_imp)
for (f in 1:length(factors)){
  f1 <- factors[f]
  model <- allModels[[f1]]
  
  # Get features for model fitting
  features <- colnames(model$trainingData)[-10001]
  
  # Make predictions
  predictedScore_factors[,f] <- predict(model, X_EMIF_imp[,features], type = "prob")$No
  
}
predictedScore_factors <- as.data.frame(predictedScore_factors)

# Predict age
X_EMIF_imp_B <- (2^X_EMIF_imp)/(1+ 2^X_EMIF_imp)
predictedAge <- agep(t(X_EMIF_imp_B), 
                     method='all')

predictedScore_factors$Age <- predictedAge$skinblood.skinblood.age

# Save predicted factors
save(predictedScore_factors, file = "EMIF/predictedScore_factors_EMIF.RData")
