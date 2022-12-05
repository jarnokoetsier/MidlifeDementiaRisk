# Clear workspace and console
rm(list = ls())
cat("\014") 

library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
source("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA/MachineLearning/FUN_MachineLearning.R")

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Load Pheno data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Load S-scored data
files <- list.files('X_S')
for (f in files){
  load(paste0("X_S/",f))
}
# Load var data
files <- list.files('X_var')
for (f in files){
  load(paste0("X_var/",f))
}

# Load varM data
files <- list.files('X_varM')
for (f in files){
  load(paste0("X_varM/",f))
}


#*****************************************************************************#
# Model training
#*****************************************************************************#

#=============================================================================#
# FILL IN
#=============================================================================#

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "var"

# Prepare data (M-values)
X_train = log2(X_CAIDE1_var/(1-X_CAIDE1_var))
Y_train = Y_CAIDE1$CAIDE

# Test if samples are in correct order
all(colnames(X_train) == Y_CAIDE1$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5

#=============================================================================#

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           number = nfold, 
                           repeats = nrep, 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = regressionSummary)

# Set grid for lambda (2.5 for CAIDE)
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "spls"

# Register cores for parallel computing
#detectCores()
#nCores <- 3
#cl <- makeCluster(nCores)
#registerDoParallel(cl)

# Actual training
set.seed(123)
fit <- train(x = t(X_train),
             y = Y_train,
             metric= performance_metric,
             method = MLmethod,
             trControl = fitControl,
             maximize = FALSE)

# Stop clusters
#stopCluster(cl)

# Get results
trainResults <- fit$results

# Get optimal lambda and alpha
optAlpha <- fit$bestTune$alpha
optLambda <- fit$bestTune$lambda

# Get coefficient value during the repeated CV
coefs <- matrix(NA, ncol(t(X_train))+1, nfold*nrep)
perf <- rep(NA, nfold*nrep)
folds <- fit$control$index
pred_CV <- NULL
obs_CV <- NULL
fold_CV <- NULL
count = 0
for (r in 1:nrep){
  for (f in 1:nfold){
    count = count + 1
    en_model_cv <- glmnet(x = t(X_train)[folds[[count]],], 
                          y = Y_train[folds[[count]]], 
                          family = "gaussian",
                          alpha = optAlpha, 
                          lambda = optLambda,
                          standardize = TRUE)
    
    # Get coefficients
    coefs[,count] <- as.matrix(coef(en_model_cv))
    
    # Get prediction
    pred <- predict(en_model_cv, t(X_train)[-folds[[count]],])[,1]
    
    # Get performance
    perf[count] <- RMSE(pred = pred, obs = Y_train[-folds[[count]]])
    
    pred_CV <- c(pred_CV,pred)
    obs_CV <- c(obs_CV, Y_train[-folds[[count]]])
    fold_CV <- c(fold_CV, rep(count,length(pred)))
  }
}

ObsPred_CV <- data.frame(Predicted = pred_CV,
                         Observed = obs_CV,
                         Fold = fold_CV)

# Get final model
finalModel <- glmnet(x = t(X_train), 
                     y = Y_train, 
                     family = "gaussian",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)

# Save output
save(trainResults, optLambda, optAlpha, perf, ObsPred_CV, coefs, finalModel,
     file = paste0("CV_", Score, "_", FeatureSelection,".RData"))

