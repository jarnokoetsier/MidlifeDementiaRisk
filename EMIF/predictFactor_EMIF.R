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

save(predictedScore_factors, file = "EMIF/predictedScore_factors_EMIF.RData")

metaData_EMIF <- as.data.frame(metaData_EMIF)
rownames(metaData_EMIF) <- metaData_EMIF$X
samples <- intersect(metaData_EMIF$X, rownames(predictedScore_factors))

metaData_fil <- metaData_EMIF[samples,]
predictedScore_factors_fil <- predictedScore_factors[samples,]

plot(metaData_fil$Hypertension, predictedScore_factors_fil$SysBP)
cor.test(metaData_fil$Hypertension, predictedScore_factors_fil$SysBP, method = "spearman")

plot(metaData_fil$Obesity, predictedScore_factors_fil$BMI)
cor.test(metaData_fil$Obesity, predictedScore_factors_fil$BMI, method = "spearman")

plot(metaData_fil$Alcohol, predictedScore_factors_fil$Alcohol)
cor.test(metaData_fil$Alcohol, predictedScore_factors_fil$Alcohol, method = "spearman")

plot(metaData_fil$CardiovascularDis, predictedScore_factors_fil$HeartDisease)
cor.test(metaData_fil$CardiovascularDis, predictedScore_factors_fil$HeartDisease, method = "spearman")

plot(metaData_fil$Depression, predictedScore_factors_fil$Depression)
cor.test(metaData_fil$Depression, predictedScore_factors_fil$Depression, method = "spearman")

plot(metaData_fil$Eduy_IMP,  predictedScore_factors_fil$Education)
cor.test(metaData_fil$Eduy_IMP, predictedScore_factors_fil$Education, method = "spearman")

plot(metaData_fil$Eduy,  predictedScore_factors_fil$Education)
cor.test(metaData_fil$Eduy, predictedScore_factors_fil$Education, method = "spearman")


plot(metaData_fil$SmokingAmount,  predictedScore_factors_fil$Smoking)
cor.test(as.numeric(metaData_fil$SmokingAmount),  predictedScore_factors_fil$Smoking, method = "spearman")

plot(metaData_fil$CholesterolloweringDrugs,  predictedScore_factors_fil$HDL)
cor.test(metaData_fil$CholesterolloweringDrugs,  predictedScore_factors_fil$HDL, method = "spearman")

library(pROC)
test <- roc(metaData_fil$Depression, predictedScore_factors_fil$Depression)
