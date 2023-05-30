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
load("~/EMIF/lumi_dpval_EMIF.RData")
load("~/EMIF/X_EMIF.RData")
load("~/EMIF/metaData_EMIF.RData")
#load("~/allModels.RData")

# optimize imputation
X_EMIF_M <- log2(X_EMIF/(1-X_EMIF))
test <- lumi_dpval[rownames(X_EMIF_M), colnames(X_EMIF_M)]
X_EMIF_mis <- X_EMIF_M
X_EMIF_mis[test > 0.1] <- NA

# remove random probes
set.seed(123)
removeProbes <- matrix(NA, nrow = 100, ncol = 2)
for (i in 1:100){
  removeProbes[i,1] <- sample(1:nrow(X_EMIF_mis),1)
  removeProbes[i,2] <- sample(1:ncol(X_EMIF_mis),1)
}
removeProbes <- removeProbes[!duplicated(removeProbes),]

# Collect true values of the probes
X_EMIF_copy <- X_EMIF_mis
values <- rep(NA, 100)
for (j in 1:100){
  values[j] <- X_EMIF_mis[removeProbes[j,1], removeProbes[j,2]]
  X_EMIF_copy[removeProbes[j,1], removeProbes[j,2]] <- NA
}


pred <- matrix(NA, nrow = 100,ncol = 10)
nPCs_CV <- 20
for (p in 15:nPCs_CV){
  
  # Impute missing values
  set.seed(123)
  X_EMIF_imp_CV <- imputePCA(t(X_EMIF_copy),ncp = p)$completeObs
  
  # Get predicted values of randomly removed probes
  for (j in 1:100){
    pred[j,p] <- X_EMIF_imp_CV[removeProbes[j,2], removeProbes[j,1]]
  }
  # Save predicted value
  save(pred, file = "EMIF/pred_imp2.RData")
}


load("EMIF/pred_imp.RData")

# Evaluate which number of PCs is the best
test <- apply(pred,2,function(x) RMSE(obs = values,pred = x))
optPC <- which.min(test) # 15 PCs

# Perform imputation using optimal number of PCs
X_EMIF_imp <- imputePCA(t(X_EMIF_mis),ncp = optPC)$completeObs
save(X_EMIF_imp, file = "EMIF/X_EMIF_imp.RData")









