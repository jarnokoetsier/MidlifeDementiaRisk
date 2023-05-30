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

# Set direcotries
DataDir <- "Data/"
CorDir <- "CAIDE1_Cor/"
OutputDir <- "PreProcessing/"

# Load data
load(paste0(DataDir,"CVindex_CAIDE1.RData"))

load(paste0(CorDir,"X_CAIDE1_CorCV.RData"))
load(paste0(CorDir,"X_CAIDE1_Cor.RData"))
load(paste0(CorDir,"X_nonTest_Cor.RData"))
load(paste0(CorDir,"X_test_Cor.RData"))

load(paste0(DataDir,"Y_CAIDE1.RData"))
load(paste0(DataDir,"Y_nonTest.RData"))
load(paste0(DataDir,"Y_test.RData"))
load(paste0(DataDir,"cellType.RData"))

# Load machine learning functions
source("FUN_MachineLearning.R")

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

# Machine learning method
MLmethod = "glmnet"

trainResults <- list()
for (i in 1:length(CVindex)){
  
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
  
  
 trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "trainResults_LowRisk_CAIDE1_Cor_EN.RData")

# Get performance
perf <- matrix(NA, nrow = 1000, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$ROC
}
optPar <- which.max(rowMeans(perf))

# Get optimal parameters
optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]


# Load data for final model training
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

# Train final model
fitControl <- trainControl(method = "none", classProbs = TRUE)
finalModel <- train(x = t(X_train),
                    y = Y_train$Class,
                    metric= "ROC",
                    method = "glmnet",
                    tuneGrid = data.frame(
                      .alpha = optAlpha,
                      .lambda = optLambda
                    ),
                    trControl = fitControl,
                    maximize = TRUE)


# Save data
save(trainResults, optLambda, optAlpha, perf, finalModel,
     file = paste0("CV_CAIDE1_Cor_LowRisk_EN.RData"))



X_test <- log2(X_test_Cor/(1-X_test_Cor))
test <- predict(finalModel, t(X_test), type = "prob")

plot(ifelse(Y_test$CAIDE < 4, 2, 1),test[,1])


################################################################################
#
# Model training CAIDE1 (EN) -HighRisk
#
################################################################################

# Class
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Intermediate_Low"),
                        levels = c("Intermediate_Low","High"))

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Machine learning method
MLmethod = "glmnet"


trainResults <- list()
for(i in 1:length(CVindex)){
  
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
  
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "trainResults_HighRisk_CAIDE1_Cor_EN.RData")


# Get performance
perf <- matrix(NA, nrow = 1000, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$ROC
}
optPar <- which.max(rowMeans(perf))

# Get optimal parameters
optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

# Load data for final model training
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

# Train final model
fitControl <- trainControl(method = "none", classProbs = TRUE)
finalModel <- train(x = t(X_train),
                    y = Y_train$Class,
                    metric= "ROC",
                    method = "glmnet",
                    tuneGrid = data.frame(
                      .alpha = optAlpha,
                      .lambda = optLambda
                    ),
                    trControl = fitControl,
                    maximize = TRUE)


# Save data
save(trainResults, optLambda, optAlpha, perf, finalModel,
     file = paste0("CV_CAIDE1_Cor_HighRisk_EN.RData"))


################################################################################
#
# Model training CAIDE1 (sPLSDA) -LowRisk
#
################################################################################

Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE < 4, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Number of component (K)
K_CV <- 1:20

# Thresholding parameter (eta)
eta_CV <- seq(0.1,0.9,length.out = 20)

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = 0.5

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

# Machine learning method
MLmethod = "spls"

trainResults <- list()
for (i in 1:length(CVindex)){
  
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
  
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "trainResults_LowRisk_CAIDE1_Cor_sPLSDA.RData")

# Get performance
perf <- matrix(NA, nrow = 400, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$ROC
}
optPar <- which.max(rowMeans(perf))

# Get optimal parameters
optK <- trainResults[[1]]$K[optPar]
optEta <- trainResults[[1]]$eta[optPar]
optKappa <- trainResults[[1]]$kappa[optPar]


# Load data for final model training
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

# Train final model
fitControl <- trainControl(method = "none", classProbs = TRUE)
finalModel <- train(x = t(X_train),
                    y = Y_train$Class,
                    metric= "ROC",
                    method = "spls",
                    tuneGrid = data.frame(
                      .K = optK,
                      .eta = optEta,
                      .kappa = optKappa
                    ),
                    trControl = fitControl,
                    maximize = TRUE)


# Save data
save(trainResults, trainResults, optK, optEta, optKappa, perf, finalModel,
     file = paste0("CAIDE1_Cor/CV_CAIDE1_Cor_LowRisk_sPLSDA.RData"))



X_test <- log2(X_test_Cor/(1-X_test_Cor))
test <- predict(finalModel, t(X_test), type = "prob")

plot(ifelse(Y_test$CAIDE < 4, 2, 1),test[,1])


################################################################################
#
# Model training CAIDE1 (sPLSDA) -HighRisk
#
################################################################################


# Class
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Intermediate_Low"),
                        levels = c("Intermediate_Low","High"))

# Number of component (K)
K_CV <- 1:20

# Thresholding parameter (eta)
eta_CV <- seq(0.1,0.9,length.out = 20)

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = 0.5

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

# Machine learning method
MLmethod = "spls"

trainResults <- list()
for(i in 1:length(CVindex)){
  
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
  
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "CAIDE1_Cor/trainResults_HighRisk_CAIDE1_Cor_sPLSDA.RData")


perf <- matrix(NA, nrow = 400, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$ROC
}
optPar <- which.max(rowMeans(perf))

optK <- trainResults[[1]]$K[optPar]
optEta <- trainResults[[1]]$eta[optPar]
optKappa <- trainResults[[1]]$kappa[optPar]

# Load data for final model training
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

# Train final model
fitControl <- trainControl(method = "none", classProbs = TRUE)
finalModel <- train(x = t(X_train),
                    y = Y_train$Class,
                    metric= "ROC",
                    method = "spls",
                    tuneGrid = data.frame(
                      .K = optK,
                      .eta = optEta,
                      .kappa = optKappa
                    ),
                    trControl = fitControl,
                    maximize = TRUE)


# Save data
save(trainResults, optK, optEta, optKappa, perf, finalModel,
     file = paste0("CAIDE1_Cor/CV_CAIDE1_Cor_HighRisk_sPLSDA.RData"))

################################################################################
#
# Model training CAIDE1 (RF) -LowRisk
#
################################################################################

library(e1071)
library(ranger)
library(dplyr)

Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE < 4, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Number of randomly selected predictors
mtry_CV <- c(100, 500,1000,1500,2000,3000,4000)

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = c(3,5,10,15,20)

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
MLmethod = "ranger"

trainResults <- list()
for (i in 1:length(CVindex)){
  
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
  
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "CAIDE1_Cor/trainResults_LowRisk_CAIDE1_Cor_RF.RData")

# Get performance
perf <- matrix(NA, nrow = 35, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$ROC
}
optPar <- which.max(rowMeans(perf))

# Get optimal parameters
opt_mtry <- trainResults[[1]]$mtry[optPar]
opt_splitrule <- trainResults[[1]]$splitrule[optPar]
opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]

# Load data for final model training
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

# Train final model
fitControl <- trainControl(method = "none", classProbs = TRUE)
set.seed(123)
finalModel <- train(x = t(X_train),
                    y = Y_train$Class,
                    metric= "ROC",
                    method = "ranger",
                    tuneGrid = data.frame(
                      .mtry = opt_mtry, 
                      .splitrule = opt_splitrule, 
                      .min.node.size = opt_min.node.size
                    ),
                    trControl = fitControl,
                    maximize = TRUE)


# Save data
save(trainResults, trainResults, opt_mtry, opt_splitrule, opt_min.node.size, perf, finalModel,
     file = paste0("CAIDE1_Cor/CV_CAIDE1_Cor_LowRisk_RF.RData"))


################################################################################
#
# Model training CAIDE1 (RF) -HighRisk
#
################################################################################

library(e1071)
library(ranger)
library(dplyr)


Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))


# Number of randomly selected predictors
mtry_CV <- c(100, 500,1000,1500,2000,3000,4000)

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = c(3,5,10,15,20)

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
MLmethod = "ranger"

trainResults <- list()
for (i in 1:length(CVindex)){
  
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
  
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "trainResults_HighRisk_CAIDE1_Cor_RF.RData")

# Get performance
perf <- matrix(NA, nrow = 35, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$ROC
}
optPar <- which.max(rowMeans(perf))

# Get optimal parameters
opt_mtry <- trainResults[[1]]$mtry[optPar]
opt_splitrule <- trainResults[[1]]$splitrule[optPar]
opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]

# Load data for final model training
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
all(colnames(X_train) == Y_train$Basename)

# Train final model
fitControl <- trainControl(method = "none", classProbs = TRUE)
finalModel <- train(x = t(X_train),
                    y = Y_train$Class,
                    metric= "ROC",
                    method = "ranger",
                    tuneGrid = data.frame(
                      .mtry = opt_mtry,
                      .splitrule = opt_splitrule,
                      .min.node.size = opt_min.node.size
                    ),
                    trControl = fitControl,
                    maximize = TRUE)


# Save data
save(trainResults, trainResults, opt_mtry, opt_splitrule, opt_min.node.size, perf, finalModel,
     file = paste0("CAIDE1_Cor/CV_CAIDE1_Cor_HighRisk_RF.RData"))

