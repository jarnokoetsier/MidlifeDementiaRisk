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
CorDir <- "LIBRA_Cor/"

# Load data
load(paste0(DataDir,"CVindex_LIBRA.RData"))

load(paste0(CorDir,"X_LIBRA_CorLCV.RData"))
load(paste0(CorDir,"X_LIBRA_CorL.RData"))
load(paste0(CorDir,"X_nonTest_CorL.RData"))
load(paste0(CorDir,"X_test_CorL.RData"))

load(paste0(DataDir,"Y_LIBRA.RData"))
load(paste0(DataDir,"Y_nonTest.RData"))
load(paste0(DataDir,"Y_test.RData"))
load(paste0(DataDir,"cellType.RData"))

# Load machine learning functions
source("FUN_MachineLearning.R")

# Prepare data
X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))
Y_train = Y_LIBRA

# Test if samples are in correct order
all(colnames(X_train) == Y_LIBRA$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5


################################################################################
#
# ElasticNet
#
################################################################################

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "glmnet"

trainResults <- list()
for (i in 1:length(CVindex)){
  
  # Select samples from specific fold
  index <- list(CVindex[[i]])
  X_CV <- X_train[,index[[1]]]
  
  # Calculate correlations (using X_CV)
  factors <- Y_train[index[[1]],14:24]
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
  for (p in 1:ncol(correlations_CV)){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1100))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < ncol(correlations_CV)){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(method = "repeatedcv", 
                             search = "grid", 
                             savePredictions = FALSE,
                             summaryFunction = regressionSummary,
                             index = index)
  
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(X_train[finalProbes,]),
               y = Y_train$LIBRA,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = FALSE)
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "LIBRA_Cor/trainResults_LIBRA_Cor_EN.RData")


#==============================================================================#
# Train final model
#==============================================================================#

# Get performance
perf <- matrix(NA, nrow = 1000, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
# Optimal performance (minimum mean RMSE in CV)
optPar <- which.min(rowMeans(perf))

# Get optimal alpha and lambda
optAlpha <- trainResults[[1]]$alpha[optPar]
optLambda <- trainResults[[1]]$lambda[optPar]

# Prepare data
X_train = log2(X_LIBRA_CorL/(1-X_LIBRA_CorL))
Y_train = Y_LIBRA
all(colnames(X_train) == Y_train$Basename)

# Actual training
set.seed(123)
fitControl <- fitControl <- trainControl(method = "none")
finalModel <- train(x = t(X_train),
                    y = Y_train$LIBRA,
                    metric= "RMSE",
                    method = "glmnet",
                    tuneGrid = data.frame(
                      .alpha = optAlpha,
                      .lambda = optLambda
                    ),
                    trControl = fitControl,
                    maximize = FALSE)


save(trainResults, optLambda, optAlpha, perf, finalModel,
     file = "LIBRA_Cor/CV_LIBRA_Cor_EN.RData")



################################################################################
#
# sPLS
#
################################################################################

# Number of component (K)
K_CV <- 1:20

# Thresholding parameter (eta)
eta_CV <- seq(0.1,0.9,length.out = 20)

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = 0.5

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "spls"

trainResults <- list()
for (i in 1:length(CVindex)){
  
  # Select samples from specific fold
  index <- list(CVindex[[i]])
  X_CV <- X_train[,index[[1]]]
  
  # Calculate correlations (using X_CV)
  factors <- Y_train[index[[1]],14:24]
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
  for (p in 1:ncol(correlations_CV)){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1100))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < ncol(correlations_CV)){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(method = "repeatedcv", 
                             search = "grid", 
                             savePredictions = FALSE,
                             summaryFunction = regressionSummary,
                             index = index)
  
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(X_train[finalProbes,]),
               y = Y_train$LIBRA,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = FALSE)
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "LIBRA_Cor/trainResults_LIBRA_Cor_sPLS.RData")

#==============================================================================#
# Train final model
#==============================================================================#

perf <- matrix(NA, nrow = 400, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
optPar <- which.min(rowMeans(perf))

optK <- trainResults[[1]]$K[optPar]
optEta <- trainResults[[1]]$eta[optPar]
optKappa <- trainResults[[1]]$kappa[optPar]



X_train = log2(X_LIBRA_CorL/(1-X_LIBRA_CorL))
Y_train = Y_LIBRA
all(colnames(X_train) == Y_train$Basename)

# Actual training
set.seed(123)
fitControl <- fitControl <- trainControl(method = "none")
finalModel <- train(x = t(X_train),
                    y = Y_train$LIBRA,
                    metric= "RMSE",
                    method = "spls",
                    tuneGrid = data.frame(
                      .K = optK,
                      .eta = optEta,
                      .kappa = optKappa
                    ),
                    trControl = fitControl,
                    maximize = FALSE)


save(trainResults, optK, optEta, optKappa, perf, finalModel,
     file ="LIBRA_Cor/CV_LIBRA_Cor_sPLS.RData")


load("LIBRA/X_test_CorL.RData")

X_test = log2(X_test_CorL/(1-X_test_CorL))
test <- predict(finalModel, t(X_test))
p <- ggplot(data.frame(pred = test,
                  obs = Y_test$LIBRA)) +
  geom_point(aes(x = obs, y = pred))

ggsave(p, file = "test.png")
plot(test, Y_test$LIBRA)

RMSE(pred = test, obs = Y_test$LIBRA)

################################################################################
#
# Random Forest
#
################################################################################
library(e1071)
library(ranger)
library(dplyr)

# Number of randomly selected predictors
mtry_CV <- c(1000, 2000,3000,4000,5000,6000,7000,8000)

# split rule
splitrule_CV <- "variance"

# minimal node size
min.node.size_CV = c(10,20,30,40,50,60)

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "ranger"

trainResults <- list()
for (i in 1:length(CVindex)){
  
  # Select samples from specific fold
  index <- list(CVindex[[i]])
  X_CV <- X_train[,index[[1]]]
  
  # Calculate correlations (using X_CV)
  factors <- Y_train[index[[1]],14:24]
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
  for (p in 1:ncol(correlations_CV)){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1100))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < ncol(correlations_CV)){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(method = "repeatedcv", 
                             search = "grid", 
                             savePredictions = FALSE,
                             summaryFunction = regressionSummary,
                             index = index)
  
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(X_train[finalProbes,]),
               y = Y_train$LIBRA,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = FALSE)
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "LIBRA_Cor/trainResults_LIBRA_Cor_RF.RData")

#==============================================================================#
# Train final model
#==============================================================================#

perf <- matrix(NA, nrow = 48, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
optPar <- which.min(rowMeans(perf))

opt_mtry <- trainResults[[1]]$mtry[optPar]
opt_splitrule <- trainResults[[1]]$splitrule[optPar]
opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]



X_train = log2(X_LIBRA_CorL/(1-X_LIBRA_CorL))
Y_train = Y_LIBRA
all(colnames(X_train) == Y_train$Basename)

# Actual training
fitControl <- fitControl <- trainControl(method = "none")
set.seed(123)
finalModel <- train(x = t(X_train),
                    y = Y_train$LIBRA,
                    metric= "RMSE",
                    method = "ranger",
                    tuneGrid = data.frame(
                      .mtry = opt_mtry,
                      .splitrule = opt_splitrule,
                      .min.node.size = opt_min.node.size
                    ),
                    trControl = fitControl,
                    maximize = FALSE)


save(trainResults, opt_mtry, opt_splitrule, opt_min.node.size, perf, finalModel,
     file ="LIBRA_Cor/CV_LIBRA_Cor_RF.RData")



X_test = log2(X_test_CorL/(1-X_test_CorL))
test <- predict(finalModel, t(X_test))
p <- ggplot(data.frame(pred = test,
                       obs = Y_test$LIBRA)) +
  geom_point(aes(x = obs, y = pred))

ggsave(p, file = "test.png")
plot(test, Y_test$LIBRA)

RMSE(pred = test, obs = Y_test$CAIDE)


