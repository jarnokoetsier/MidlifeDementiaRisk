# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Load machine learning functions
source("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA/MachineLearning/FUN_MachineLearning.R")

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Load CV index
load("CVindex_CAIDE1.RData")
#=============================================================================#
# FILL IN
#=============================================================================#

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "var"

# Load data
files <- list.files(paste0("X_", FeatureSelection))
for (f in files){
  load(paste0("X_", FeatureSelection, "/", f))
}

# Prepare data
X_train = log2(X_CAIDE1_var/(1-X_CAIDE1_var))
Y_train = Y_CAIDE1$CAIDE

# Test if samples are in correct order
all(colnames(X_train) == Y_CAIDE1$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5

# Performance metric
performance_metric = "ROC"
output <- list()
#=============================================================================#



###############################################################################

# ElasticNet

###############################################################################


#*****************************************************************************#
# Model training (Low vs intermediate/high)
#*****************************************************************************#
Y_train <- factor(ifelse(Y_CAIDE1$CAIDE < 4, "Low", "Intermediate_High"),
                  levels = c("Intermediate_High","Low"))

# Settings for repeated cross-validation
fitControl <- trainControl(search = "grid", 
                           savePredictions = FALSE,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           index = CVindex)

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Machine learning method
MLmethod = "glmnet"

# Register cores for parallel computing
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Actual training
set.seed(123)
fit <- train(x = t(X_train),
             y = Y_train,
             method = MLmethod,
             family = "binomial",
             standardize = TRUE,
             metric= performance_metric,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Stop clusters
stopCluster(cl)

# Get results
trainResults <- fit$results

# Get optimal lambda and alpha
optAlpha <- fit$bestTune$alpha
optLambda <- fit$bestTune$lambda

# Get coefficients, prediction, and performance during the repeated CV
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
                          family = "binomial",
                          alpha = optAlpha, 
                          lambda = optLambda,
                          standardize = TRUE)
    
    # Get coefficients
    coefs[,count] <- as.matrix(coef(en_model_cv))
    
    # Get prediction
    pred <- predict(en_model_cv, t(X_train)[-folds[[count]],], type = "response")[,1]
    
    pred_CV <- c(pred_CV,pred)
    obs_CV <- c(obs_CV, Y_train[-folds[[count]]])
    fold_CV <- c(fold_CV, rep(count,length(pred)))
  }
}

# Save observed and predicted in a dataframe
ObsPred_CV <- data.frame(Predicted = pred_CV,
                         Observed = obs_CV,
                         Fold = fold_CV)

# Get final model
finalModel <- glmnet(x = t(X_train), 
                     y = Y_train, 
                     family = "binomial",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)

# Save output
save(trainResults, optLambda, optAlpha, ObsPred_CV, coefs, finalModel,
     file = paste0("CV_", Score, "_", FeatureSelection,"_LowRisk.RData"))

#*****************************************************************************#
# Model training (High vs intermediate/low)
#*****************************************************************************#
Y_train <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Intermediate_Low"),
                  levels = c("Intermediate_Low", "High"))

# Settings for repeated cross-validation
fitControl <- trainControl(search = "grid", 
                           savePredictions = FALSE,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           index = CVindex)

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Machine learning method
MLmethod = "glmnet"

# Register cores for parallel computing
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Actual training
set.seed(123)
fit <- train(x = t(X_train),
             y = Y_train,
             method = MLmethod,
             family = "binomial",
             standardize = TRUE,
             metric= performance_metric,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Stop clusters
stopCluster(cl)

# Get results
trainResults <- fit$results

# Get optimal lambda and alpha
optAlpha <- fit$bestTune$alpha
optLambda <- fit$bestTune$lambda

# Get coefficients, prediction, and performance during the repeated CV
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
                          family = "binomial",
                          alpha = optAlpha, 
                          lambda = optLambda,
                          standardize = TRUE)
    
    # Get coefficients
    coefs[,count] <- as.matrix(coef(en_model_cv))
    
    # Get prediction
    pred <- predict(en_model_cv, t(X_train)[-folds[[count]],], type = "response")[,1]
    
    pred_CV <- c(pred_CV,pred)
    obs_CV <- c(obs_CV, Y_train[-folds[[count]]])
    fold_CV <- c(fold_CV, rep(count,length(pred)))
  }
}

# Save observed and predicted in a dataframe
ObsPred_CV <- data.frame(Predicted = pred_CV,
                         Observed = obs_CV,
                         Fold = fold_CV)

# Get final model
finalModel <- glmnet(x = t(X_train), 
                     y = Y_train, 
                     family = "binomial",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)

# Save output
save(trainResults, optLambda, optAlpha, ObsPred_CV, coefs, finalModel,
     file = paste0("CV_", Score, "_", FeatureSelection,"_HighRisk.RData"))











