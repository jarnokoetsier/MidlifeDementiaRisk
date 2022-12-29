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
FeatureSelection = "Cor"

# Load data
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/", f))
}

# Prepare data
X_train = log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))
Y_train = Y_CAIDE1

# Test if samples are in correct order
all(colnames(X_train) == Y_CAIDE1$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5

# Performance metric
performance_metric = "ROC"
output <- list()
#=============================================================================#




################################################################################
#
# Model training CAIDE1 (EN) -LowRisk
#
################################################################################

Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE < 4, "Low", "Intermediate_High"),
                  levels = c("Intermediate_High","Low"))

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"

# Register cores for parallel computing
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

trainResults <- foreach::foreach(i = 1:length(CVindex)) %dopar% {
  
  # Select samples from specific fold
  index <- list(CVindex[[i]])
  X_CV <- X_train[,index[[1]]]
  
  # Calculate correlations (using X_CV)
  factors <- Y_train[index[[1]],14:20]
  correlations_CV <- matrix(NA, nrow = nrow(X_CV), ncol = ncol(factors))
  for (f in 1:ncol(factors)) {
    correlations_CV[,f] <- apply(X_CV, 1, 
                                 function(x){cor(x, 
                                                 factors[,f], 
                                                 method = "spearman")})
  }
  rownames(correlations_CV) <- rownames(X_CV)
  colnames(correlations_CV) <- colnames(factors)
  
  # Select top correlated features for each factor
  probes <- list()
  for (p in 1:7){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < 7){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = FALSE,
                             summaryFunction = twoClassSummary,
                             classProbs = TRUE,
                             index = index)
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(X_train[finalProbes,]),
               y = Y_train$Class,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = TRUE)

  
  return(fit$results)
  
}

# Stop clusters
stopCluster(cl)
#==============================================================================#
# Train final model
#==============================================================================#

perf <- matrix(NA, nrow = 1000, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
optPar <- which.min(rowMeans(perf))

optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

load("X_CAIDE1_Cor.RData")
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

finalModel <- glmnet(x = t(X_train), 
                     y = Y_train$Class, 
                     family = "binomial",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)

Score = "CAIDE1"
FeatureSelection = "Cor"
save(trainResults, optLambda, optAlpha, perf, finalModel,
     file = paste0("CV_", Score, "_", FeatureSelection,"_LowRisk.RData"))


################################################################################
#
# Model training CAIDE1 (EN) -HighRisk
#
################################################################################


Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Intermediate_Low"),
                        levels = c("Intermediate_Low","High"))

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"

# Register cores for parallel computing
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

trainResults <- foreach::foreach(i = 1:length(CVindex)) %dopar% {
  
  # Select samples from specific fold
  index <- list(CVindex[[i]])
  X_CV <- X_train[,index[[1]]]
  
  # Calculate correlations (using X_CV)
  factors <- Y_train[index[[1]],14:20]
  correlations_CV <- matrix(NA, nrow = nrow(X_CV), ncol = ncol(factors))
  for (f in 1:ncol(factors)) {
    correlations_CV[,f] <- apply(X_CV, 1, 
                                 function(x){cor(x, 
                                                 factors[,f], 
                                                 method = "spearman")})
  }
  rownames(correlations_CV) <- rownames(X_CV)
  colnames(correlations_CV) <- colnames(factors)
  
  # Select top correlated features for each factor
  probes <- list()
  for (p in 1:7){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < 7){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = FALSE,
                             summaryFunction = twoClassSummary,
                             classProbs = TRUE,
                             index = index)
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(X_train[finalProbes,]),
               y = Y_train$Class,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = TRUE)
  
  
  return(fit$results)
  
}

# Stop clusters
stopCluster(cl)

#==============================================================================#
# Train final model
#==============================================================================#

perf <- matrix(NA, nrow = 1000, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
optPar <- which.min(rowMeans(perf))

optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

load("X_CAIDE1_Cor.RData")
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

finalModel <- glmnet(x = t(X_train), 
                     y = Y_train$Class, 
                     family = "binomial",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)

Score = "CAIDE1"
FeatureSelection = "Cor"
save(trainResults, optLambda, optAlpha, perf, finalModel,
     file = paste0("CV_", Score, "_", FeatureSelection,"_HighRisk.RData"))









################################################################################
#
# Model training CAIDE1 (sPLSDA) -LowRisk
#
################################################################################

Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE < 4, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "spls"

# Register cores for parallel computing
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

trainResults <- foreach::foreach(i = 1:length(CVindex)) %dopar% {
  
  # Select samples from specific fold
  index <- list(CVindex[[i]])
  X_CV <- X_train[,index[[1]]]
  
  # Calculate correlations (using X_CV)
  factors <- Y_train[index[[1]],14:20]
  correlations_CV <- matrix(NA, nrow = nrow(X_CV), ncol = ncol(factors))
  for (f in 1:ncol(factors)) {
    correlations_CV[,f] <- apply(X_CV, 1, 
                                 function(x){cor(x, 
                                                 factors[,f], 
                                                 method = "spearman")})
  }
  rownames(correlations_CV) <- rownames(X_CV)
  colnames(correlations_CV) <- colnames(factors)
  
  # Select top correlated features for each factor
  probes <- list()
  for (p in 1:7){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < 7){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = FALSE,
                             summaryFunction = twoClassSummary,
                             classProbs = TRUE,
                             index = index)
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(X_train[finalProbes,]),
               y = Y_train$Class,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = TRUE)
  
  
  return(fit$results)
  
}

# Stop clusters
stopCluster(cl)
#==============================================================================#
# Train final model
#==============================================================================#

perf <- matrix(NA, nrow = 1000, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
optPar <- which.min(rowMeans(perf))

optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

load("X_CAIDE1_Cor.RData")
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

finalModel <- glmnet(x = t(X_train), 
                     y = Y_train$Class, 
                     family = "binomial",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)

Score = "CAIDE1"
FeatureSelection = "Cor"
save(trainResults, optLambda, optAlpha, perf, finalModel,
     file = paste0("CV_", Score, "_", FeatureSelection,"_LowRisk.RData"))


################################################################################
#
# Model training CAIDE1 (EN) -HighRisk
#
################################################################################


Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Intermediate_Low"),
                        levels = c("Intermediate_Low","High"))

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"

# Register cores for parallel computing
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

trainResults <- foreach::foreach(i = 1:length(CVindex)) %dopar% {
  
  # Select samples from specific fold
  index <- list(CVindex[[i]])
  X_CV <- X_train[,index[[1]]]
  
  # Calculate correlations (using X_CV)
  factors <- Y_train[index[[1]],14:20]
  correlations_CV <- matrix(NA, nrow = nrow(X_CV), ncol = ncol(factors))
  for (f in 1:ncol(factors)) {
    correlations_CV[,f] <- apply(X_CV, 1, 
                                 function(x){cor(x, 
                                                 factors[,f], 
                                                 method = "spearman")})
  }
  rownames(correlations_CV) <- rownames(X_CV)
  colnames(correlations_CV) <- colnames(factors)
  
  # Select top correlated features for each factor
  probes <- list()
  for (p in 1:7){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < 7){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = FALSE,
                             summaryFunction = twoClassSummary,
                             classProbs = TRUE,
                             index = index)
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(X_train[finalProbes,]),
               y = Y_train$Class,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = TRUE)
  
  
  return(fit$results)
  
}

# Stop clusters
stopCluster(cl)

#==============================================================================#
# Train final model
#==============================================================================#

perf <- matrix(NA, nrow = 1000, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
optPar <- which.min(rowMeans(perf))

optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

load("X_CAIDE1_Cor.RData")
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

finalModel <- glmnet(x = t(X_train), 
                     y = Y_train$Class, 
                     family = "binomial",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)

Score = "CAIDE1"
FeatureSelection = "Cor"
save(trainResults, optLambda, optAlpha, perf, finalModel,
     file = paste0("CV_", Score, "_", FeatureSelection,"_HighRisk.RData"))