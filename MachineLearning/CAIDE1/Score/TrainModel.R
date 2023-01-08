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
FeatureSelection = "CorKS"

# Load data
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/", f))
}

# Prepare data
X_train = log2(X_CAIDE1_CorKS/(1-X_CAIDE1_CorKS))
Y_train = Y_CAIDE1$CAIDE

# Test if samples are in correct order
all(colnames(X_train) == Y_CAIDE1$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5

# Performance metric
performance_metric = "RMSE"

#=============================================================================#


###############################################################################

# ElasticNet

###############################################################################


#*****************************************************************************#
# Model training
#*****************************************************************************#

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           search = "grid", 
                           savePredictions = FALSE,
                           summaryFunction = regressionSummary,
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
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

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

# Save observed and predicted in a dataframe
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



###############################################################################

# Sparse Partial Least Squares (sPLS)

###############################################################################

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           number = nfold, 
                           repeats = nrep, 
                           search = "grid", 
                           savePredictions = FALSE,
                           summaryFunction = regressionSummary)

# Number of component (K)
K_CV <- 1:10

# Thresholding parameter (eta)
eta_CV <- seq(0.1,0.9,length.out = 10)

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = 0.5

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

# Machine learning method
MLmethod = "spls"

# Standardize data
X_train_scaled <- (X_train - rowMeans(X_train))/(apply(X_train, 1, sd))

# Register cores for parallel computing
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Actual training
set.seed(123)
fit <- train(x = t(X_train_scaled),
             y = Y_train,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Stop clusters
stopCluster(cl)

# Variable importance
varImp(fit, scale = FALSE)


###############################################################################

# Extreme gradient boosting (tree)

###############################################################################
#library(xgboost)

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           number = nfold, 
                           repeats = nrep, 
                           search = "grid", 
                           savePredictions = FALSE,
                           summaryFunction = regressionSummary)

# Number of rounds (boosting iterations)
nrounds_CV <- 2:5

# Learning rate (eta)
eta_CV <- 0.3 # between 0 and 1 (default = 0.3)
# scale the contribution of each tree
# Smaller: robust to overfitting but slow to compute (more nrounds)

# Max tree depth
max.depth_CV <- 2:10 # Default = 6

# Gamma: loss of reduction required to make a further partition
# Larger: more conservative
gamma_CV <- 1

# Subsample percentage
subsample_CV <- 0.5

# subsample fo columns
colsample_bytree_CV <- 0.5

# Minimum number of instances needed to be in each node
min_child_weight_CV <- 5:20

# Combine into a single data frame
parameterGrid <- expand.grid(nrounds_CV,
                             max.depth_CV,
                             eta_CV,
                             gamma_CV,
                             colsample_bytree_CV,
                             min_child_weight_CV,
                             subsample_CV)

colnames(parameterGrid) <- c(".nrounds", 
                             ".max.depth", 
                             ".eta",
                             ".gamma",
                             ".colsample_bytree",
                             ".min_child_weight",
                             ".subsample")

# Machine learning method
MLmethod = "xgbTree"

# Register cores for parallel computing
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Actual training
set.seed(123)
fit <- train(x = t(X_train),
             y = Y_train,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Stop clusters
stopCluster(cl)

# Variable importance
varImp(fit, scale = FALSE)


