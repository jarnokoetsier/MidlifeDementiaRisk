library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(wateRmelon)
library(missMDA)

# Clear workspace and console
rm(list = ls())
cat("\014") 


# Load data
load("~/ADNI/ADNI_metadata_classes_640.Rdata")
load("ADNI/MetaData_ADNI.RData")
load("~/allModels.RData")
load("ADNI/X_ADNI_imp.RData")

# Predict factors
factors <- names(allModels)
predictedScore_factors <- matrix(NA, nrow = nrow(X_ADNI_imp), ncol = length(factors))
colnames(predictedScore_factors) <- factors
rownames(predictedScore_factors) <- rownames(X_ADNI_imp)
for (f in 1:length(factors)){
  f1 <- factors[f]
  model <- allModels[[f1]]
  
  # Get features for model fitting
  features <- colnames(model$trainingData)[-10001]
  
  # Make predictions
  predictedScore_factors[,f] <- predict(model, X_ADNI_imp[,features], type = "prob")$No
  
}
predictedScore_factors <- as.data.frame(predictedScore_factors)


# Predict age
X_ADNI_imp_B <- (2^X_ADNI_imp)/(1+ 2^X_ADNI_imp)
predictedAge <- agep(t(X_ADNI_imp_B), 
                     method='all')

predictedScore_factors$Age <- predictedAge$skinblood.skinblood.age

save(predictedScore_factors, file = "ADNI/predictedScore_factors_ADNI.RData")



# Accelated age
MetaData_baseline <- as.data.frame(MetaData_baseline)
rownames(MetaData_baseline) <- MetaData_baseline$Basename
MetaData_predictedScore <- MetaData_baseline[rownames(predictedScore_factors),]

plot(MetaData_predictedScore$Age, predictedScore_factors$Age)

getAge <- data.frame(obsAge = MetaData_predictedScore$Age,
                     predAge = predictedScore_factors$Age)
rownames(getAge) <- rownames(predictedScore_factors)
ageModel <- lm(obsAge ~ predAge, data = getAge)
accAge <- residuals(ageModel)

all(names(accAge) == rownames(predictedScore_factors))
predictedScore_factors$AccAge <- accAge





Y_ADNI <- left_join(MetaData_predictedScore, ADNI_model_res,
                    by = c("RID" = "RID"))

save(Y_ADNI, file = "ADNI/Y_ADNI.RData")

length(unique(MetaData_predictedScore$RID))








# Keep midlife samples only:
keepSamples <- MetaData_baseline$Basename[MetaData_baseline$Age <= 75]
keepSamples <- colnames(X_ADNI_M)


X_ADNI_fil <- X_ADNI_M[,colnames(X_ADNI_M) %in%keepSamples]
