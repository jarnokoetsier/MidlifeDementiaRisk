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
CorDir <- "CAIDE2_Cor/"

# Load data
load(paste0(DataDir,"CVindex_CAIDE2.RData"))

load(paste0(CorDir,"X_CAIDE2_Cor2CV.RData"))
load(paste0(CorDir,"X_CAIDE2_Cor2.RData"))
load(paste0(CorDir,"X_nonTest_Cor2.RData"))
load(paste0(CorDir,"X_test_Cor2.RData"))

load(paste0(DataDir,"Y_CAIDE2.RData"))
load(paste0(DataDir,"Y_nonTest.RData"))
load(paste0(DataDir,"Y_test.RData"))
load(paste0(DataDir,"cellType.RData"))

# Load machine learning functions
source("FUN_MachineLearning.R")

# Prepare data
X_train = log2(X_CAIDE2_Cor2CV/(1-X_CAIDE2_Cor2CV))
Y_train = Y_CAIDE2

# Test if samples are in correct order
all(colnames(X_train) == Y_CAIDE2$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5

################################################################################
#
# Model training CAIDE2 (EN)
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
  for (p in 1:ncol(correlations_CV)){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1650))
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
               y = Y_train$CAIDE2,
               metric= performance_metric,
               method = MLmethod,
               tuneGrid = parameterGrid,
               trControl = fitControl,
               maximize = FALSE)
  
  trainResults[[i]] <- fit$results
  
}
save(trainResults, file = "CAIDE2_Cor/trainResults_CAIDE2_Cor_EN.RData")


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
X_train = log2(X_CAIDE2_Cor2/(1-X_CAIDE2_Cor2))
Y_train = Y_CAIDE2
all(colnames(X_train) == Y_train$Basename)

# Actual training
set.seed(123)
fitControl <- fitControl <- trainControl(method = "none")
finalModel <- train(x = t(X_train),
                    y = Y_train$CAIDE2,
                    metric= "RMSE",
                    method = "glmnet",
                    tuneGrid = data.frame(
                      .alpha = optAlpha,
                      .lambda = optLambda
                    ),
                    trControl = fitControl,
                    maximize = FALSE)


save(trainResults, optLambda, optAlpha, perf, finalModel,
     file = paste0("CV_CAIDE2_Cor_EN.RData"))
