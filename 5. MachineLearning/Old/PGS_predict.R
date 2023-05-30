

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
DataDir <- "E:/Thesis/EXTEND/Methylation/"

# Load data
load(paste0(DataDir,"CVindex_CAIDE2.RData"))

# Load phenotypes
load(paste0(DataDir,"Y/Y_CAIDE2.RData"))
load(paste0(DataDir,"Y/Y_nonTest.RData"))
load(paste0(DataDir,"Y/Y_test.RData"))
load(paste0(DataDir,"cellType.RData"))

# Load PGS
load("E:/Thesis/EXTEND/df_list.RData")
dataMatrix <- df_list$`bayesr-shrink`

samples <- intersect(rownames(dataMatrix),Y_CAIDE2$ID)

rownames(Y_CAIDE2) <- Y_CAIDE2$ID
Y_train <- Y_CAIDE2[samples,]

X_train <- dataMatrix[samples,]
X_train <- t(t((t(X_train)-rowMeans(t(X_train)))/apply(t(X_train),1,sd)))


# Test if samples are in correct order
all(colnames(X_train) == Y_train$ID)

# Load machine learning functions
source("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA/MachineLearning/FUN_MachineLearning.R")


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

# Actual training
set.seed(123)
fitControl <- fitControl <- trainControl(method = "none")
finalModel <- train(x = t(X_train),
                    y = Y_train$CAIDE,
                    metric= "RMSE",
                    method = "glmnet",
                    tuneGrid = data.frame(
                      .alpha = optAlpha,
                      .lambda = optLambda
                    ),
                    trControl = fitControl,
                    maximize = FALSE)



test <- varImp(finalModel)$importance


# Linear model + RFE
feature_list <- list()
features <- colnames(X_train)
perf_all <- matrix(NA, nrow = length(CVindex), ncol = length(features))
for (f in 1:length(features)){
  performance <- rep(NA, length(CVindex))
  coeffs <- matrix(NA,nrow = length(features), ncol = length(CVindex))
  feature_list[[f]] <- features
  for (i in 1:length(CVindex)){
    
    # Select samples from specific fold
    index <- list(CVindex[[i]])
    X_CV <- as.data.frame(X_train[index[[1]],])
    Y_CV <- Y_train[index[[1]],"CAIDE"]
    
    X_val <- as.data.frame(X_train[-index[[1]],])
    Y_val <- Y_train[-index[[1]],"CAIDE"]
    
    # Fit model
    fitData <- cbind.data.frame(Y_CV, X_CV[,features])
    colnames(fitData) <- c("Y_CV", features)
    test <- lm(Y_CV ~ .,data = fitData)
    
    # Get performance
    pred <- predict(test, X_val)
    performance[i] <- RMSE(pred = pred, obs = Y_val)
    
    # Collect coefficients
    coeffs[,i] <- coef(test)[-1]
  }
  
  # Remove feature with lowest absolute coefficient
  removeFeature <- features[which.min(rowMeans(abs(coeffs)))]
  features <- setdiff(features, removeFeature)
  
  perf_all[,f] <- performance
}


