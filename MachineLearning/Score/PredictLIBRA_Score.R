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
# Model training CAIDE1 (EN)
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
X_train = log2(X_LIBRA_Cor/(1-X_LIBRA_Cor))
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

# Prepare data (M-values)
load("LIBRA/X_LIBRA_CorCV.RData")
load("LIBRA/CVindex_LIBRA.RData")

X_train = log2(X_LIBRA_CorCV/(1-X_LIBRA_CorCV))
Y_train = Y_LIBRA
all(colnames(X_train) == Y_train$Basename)

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
for (i in 4:length(CVindex)){
  
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
save(trainResults, file = "trainResults_LIBRA_Cor_sPLS.RData")

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



load("LIBRA/X_LIBRA_CorL.RData")
X_train = log2(X_LIBRA_Cor/(1-X_LIBRA_Cor))
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
     file = paste0("CV_LIBRA_Cor_sPLS.RData"))


load("LIBRA/X_test_CorL.RData")

X_test = log2(X_test_Cor/(1-X_test_Cor))
test <- predict(finalModel, t(X_test))
p <- ggplot(data.frame(pred = test,
                  obs = Y_test$LIBRA)) +
  geom_point(aes(x = obs, y = pred))

ggsave(p, file = "test.png")
plot(test, Y_test$LIBRA)

RMSE(pred = test, obs = Y_test$LIBRA)


################################################################################
#
# SVM
#
################################################################################
# Prepare data (M-values)
load("LIBRA/X_LIBRA_CorCV.RData")
load("LIBRA/CVindex_LIBRA.RData")

X_train = log2(X_LIBRA_CorCV/(1-X_LIBRA_CorCV))
Y_train = Y_LIBRA
all(colnames(X_train) == Y_train$Basename)

#sigma
#sigma_CV <- exp(seq(log(0.000001),log(0.0001),length.out = 10))
sigma_CV <- exp(seq(log(0.001),log(1),length.out = 10))

# cost
C_CV <- exp(seq(log(0.01),log(2.5),length.out = 10))

# Combine into a single data frame
parameterGrid <- expand.grid(sigma_CV, C_CV)
colnames(parameterGrid) <- c(".sigma", ".C")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "svmRadial"

trainResults <- list()
for (i in 1:3){
  
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
save(trainResults, file = "trainResults_LIBRA_Cor_SVM.RData")

#==============================================================================#
# Train final model
#==============================================================================#

perf <- matrix(NA, nrow = 100, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
optPar <- which.min(rowMeans(perf))

optSigma <- trainResults[[1]]$sigma[optPar]
optC <- trainResults[[1]]$C[optPar]




load("LIBRA/X_LIBRA_CorL.RData")
X_train = log2(X_LIBRA_Cor/(1-X_LIBRA_Cor))
Y_train = Y_LIBRA
all(colnames(X_train) == Y_train$Basename)

# Actual training
set.seed(123)
fitControl <- fitControl <- trainControl(method = "none")
finalModel <- train(x = t(X_train),
                    y = Y_train$LIBRA,
                    metric= "RMSE",
                    method = "svmRadial",
                    tuneGrid = data.frame(
                      .sigma = optSigma,
                      .C = optC
                    ),
                    trControl = fitControl,
                    maximize = FALSE)


save(trainResults, optSigma, optC, perf, finalModel,
     file = paste0("CV_LIBRA_Cor_SVM.RData"))


load("LIBRA/X_test_CorL.RData")

X_test = log2(X_test_Cor/(1-X_test_Cor))
test <- predict(finalModel, t(X_test))
p <- ggplot(data.frame(pred = test,
                       obs = Y_test$LIBRA)) +
  geom_point(aes(x = obs, y = pred))

ggsave(p, file = "test.png")
plot(test, Y_test$LIBRA)

RMSE(pred = test, obs = Y_test$LIBRA)




trainResults <- list()
for (i in 21:length(CVindex)){
  
  # Select samples from specific fold
  index <- list(CVindex[[i]])
  X_CV <- X_train[,index[[1]]]
  
  # Calculate correlations (using X_CV)
  factors <- Y_train[index[[1]],14:21]
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
  fitControl <- trainControl(method = "repeatedcv", 
                             search = "grid", 
                             savePredictions = FALSE,
                             summaryFunction = regressionSummary,
                             index = index,
                             preProc = c("center", "scale"))
  
  
  # Actual training
  set.seed(123)
  fit <- train(x = t(X_train[finalProbes,]),
               y = Y_train$CAIDE,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = FALSE)
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "trainResults_svm_cor1.RData")

#==============================================================================#
# Train final model
#==============================================================================#
load("trainResults_svm_cor1.RData")
perf <- matrix(NA, nrow = 100, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
optPar <- which.min(rowMeans(perf))

finalPar <- data.frame(
  .sigma = trainResults[[1]]$sigma[optPar],
  .C = trainResults[[1]]$C[optPar]
)


load("X_CAIDE1_Cor.RData")
X_train = log2(X_CAIDE1_Cor/(1-X_CAIDE1_Cor))
Y_train = Y_CAIDE1
all(colnames(X_train) == Y_train$Basename)

# Settings for repeated cross-validation
fitControl <- fitControl <- trainControl(method = "none")

# Actual training
set.seed(123)
finalModel <- train(x = t(X_train),
                    y = Y_train$CAIDE,
                    metric= performance_metric,
                    method = MLmethod,
                    tuneGrid = finalPar,
                    trControl = fitControl,
                    maximize = FALSE)


Score = "CAIDE1"
FeatureSelection = "Cor"
method = "svm"
save(trainResults, finalPar, perf, finalModel,
     file = paste0("CV_", Score, "_", FeatureSelection,"_",method,".RData"))

load("CV_CAIDE1_Cor_svm.RData")
load("X_test_Cor.RData")

X_test = log2(X_test_Cor/(1-X_test_Cor))
#X_test <- (X_test - means[rownames(X_test)])/(std[rownames(X_test)])
test1 <- predict(finalModel, t(X_test))
plot(test1, Y_test$CAIDE)

RMSE(pred = test1, obs = Y_test$CAIDE)