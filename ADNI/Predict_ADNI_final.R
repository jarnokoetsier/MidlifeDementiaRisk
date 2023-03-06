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
load("ADNI/X_ADNI.RData")
load("ADNI/lumi_dpval_ADNI.RData")
load("ADNI/MetaData_ADNI.RData")
load("~/allModels.RData")


# Impute missing values
X_ADNI_M <- log2(X_ADNI/(1-X_ADNI))
test <- lumi_dpval[rownames(X_ADNI_M), colnames(X_ADNI_M)]
X_ADNI_mis <- X_ADNI_M
X_ADNI_mis[test > 0.1] <- NA
X_ADNI_imp <- imputePCA(t(X_ADNI_mis),ncp = 5)$completeObs
save(X_ADNI_imp, file = "ADNI/X_ADNI_imp.RData")



# Keep midlife samples only:
keepSamples <- MetaData_baseline$Basename[MetaData_baseline$Age <= 75]
keepSamples <- colnames(X_ADNI_M)


X_ADNI_fil <- X_ADNI_M[,colnames(X_ADNI_M) %in%keepSamples]







factors <- names(allModels)
removeSamples <- NULL
for (f in 1:length(factors)){
  f1 <- factors[f]
  model <- allModels[[f1]]
  
  # Get (non-zero) coefficients from model
  coefficients <- varImp(model,scale = FALSE)$importance
  coefficients <- rownames(coefficients)[coefficients$Overall != 0]
  
  # Keep samples with low detection p-value for these features/probes
  test <- lumi_dpval[coefficients,]
  removeSamples <- c(removeSamples,colnames(test)[colSums(test>0.1) > 0])
}

X_ADNI_all <- X_ADNI_fil[,!(colnames(X_ADNI_fil) %in% removeSamples)]

factors <- names(allModels)
predictedScore <- matrix(NA, nrow = ncol(X_ADNI_all), ncol = length(factors))
colnames(predictedScore) <- factors
rownames(predictedScore) <- colnames(X_ADNI_all)
for (f in 1:length(factors)){
  f1 <- factors[f]
  model <- allModels[[f1]]
  
  # Get features for model fitting
  features <- colnames(model$trainingData)[-10001]
  
  # Make predictions
  predictedScore[,f] <- predict(model, t(X_ADNI_all[features,]), type = "prob")$No
  
}
predictedScore <- as.data.frame(predictedScore)

predictedAge <- agep(X_ADNI_all[rownames(X_ADNI_fil) %in% rownames(lumi_dpval)[rowSums(lumi_dpval>0.1) == 0],], 
                     method='all')


predictedScore$Age <- predictedAge$skinblood.skinblood.age


MetaData_baseline <- as.data.frame(MetaData_baseline)
rownames(MetaData_baseline) <- MetaData_baseline$Basename
MetaData_predictedScore <- MetaData_baseline[rownames(predictedScore),]

plot(MetaData_predictedScore$Age, predictedScore$Age)

getAge <- data.frame(obsAge = MetaData_predictedScore$Age,
                     predAge = predictedScore$Age)
rownames(getAge) <- rownames(predictedScore)
ageModel <- lm(obsAge ~ predAge, data = getAge)
accAge <- residuals(ageModel)

all(names(accAge) == rownames(predictedScore))
predictedScore$AccAge <- accAge


save(predictedScore, file = "ADNI/predictedScore.RData")

Y_ADNI <- left_join(MetaData_predictedScore, ADNI_model_res,
                     by = c("RID" = "RID"))

save(Y_ADNI, file = "ADNI/Y_ADNI.RData")

length(unique(MetaData_predictedScore$RID))
