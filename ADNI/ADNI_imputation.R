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


# optimize imputation
X_ADNI_M <- log2(X_ADNI/(1-X_ADNI))
test <- lumi_dpval[rownames(X_ADNI_M), colnames(X_ADNI_M)]
X_ADNI_mis <- X_ADNI_M
X_ADNI_mis[test > 0.1] <- NA

# remove random probes
set.seed(123)
removeProbes <- matrix(NA, nrow = 100, ncol = 2)
for (i in 1:100){
  removeProbes[i,1] <- sample(1:nrow(X_ADNI_mis),1)
  removeProbes[i,2] <- sample(1:ncol(X_ADNI_mis),1)
}
removeProbes <- removeProbes[!duplicated(removeProbes),]

# Collect true values of the probes
X_ADNI_copy <- X_ADNI_mis
values <- rep(NA, 100)
for (j in 1:100){
  values[j] <- X_ADNI_mis[removeProbes[j,1], removeProbes[j,2]]
  X_ADNI_copy[removeProbes[j,1], removeProbes[j,2]] <- NA
}


pred <- matrix(NA, nrow = 100,ncol = 20)
nPCs_CV <- 20
for (p in 1:nPCs_CV){
  
  # Impute missing values
  set.seed(123)
  X_ADNI_imp_CV <- imputePCA(t(X_ADNI_copy),ncp = p)$completeObs
  
  # Get predicted values of randomly removed probes
  for (j in 1:100){
    pred[j,p] <- X_ADNI_imp_CV[removeProbes[j,2], removeProbes[j,1]]
  }
}
save(pred, file = "ADNI/pred_imp.RData")


test <- apply(pred,2,function(x) RMSE(obs = values,pred = x))
optPC <- which.min(test) # 5 PCs
X_ADNI_imp <- imputePCA(t(X_ADNI_mis),ncp = optPC)$completeObs
save(X_ADNI_imp, file = "ADNI/X_ADNI_imp.RData")






















