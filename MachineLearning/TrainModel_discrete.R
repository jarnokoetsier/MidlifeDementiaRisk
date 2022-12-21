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
performance_metric = "Kappa"

#=============================================================================#



###############################################################################

# ElasticNet

###############################################################################

hist(Y_train)
quantile(Y_train,0.33)
quantile(Y_train,0.67)

Y_train <- ifelse(Y_train < 4, "Low", "Intermediate_High")

#*****************************************************************************#
# Model training (Low vs intermediate/high)
#*****************************************************************************#

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           number = nfold, 
                           repeats = nrep, 
                           search = "grid", 
                           savePredictions = FALSE)

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


# Get final model
finalModel <- glmnet(x = t(X_train), 
                     y = Y_train, 
                     family = "binomial",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)


pred_test <- predict(finalModel, t(X_test_var), type = "response")[,1]
response_test <- ifelse(Y_test$CAIDE < 4, 0, 1)

roc_list_test <- roc(response = response_test, 
                     predictor = pred_test)

plot(roc_list_test)

pred_test <- predict(finalModel, t(X_train), type = "response")[,1]
response_test <- ifelse(Y_CAIDE1$CAIDE < 4, 0, 1)

#*****************************************************************************#
# Model training (Low vs intermediate/high)
#*****************************************************************************#

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           number = nfold, 
                           repeats = nrep, 
                           search = "grid", 
                           savePredictions = FALSE)

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
             maximize = TRUE)

# Stop clusters
stopCluster(cl)

# Get results
trainResults <- fit$results

# Get optimal lambda and alpha
optAlpha <- fit$bestTune$alpha
optLambda <- fit$bestTune$lambda


# Get final model
finalModel <- glmnet(x = t(X_train), 
                     y = Y_train, 
                     family = "binomial",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)


pred_test <- predict(finalModel, t(X_test_cor), type = "response")[,1]
response_test <- ifelse(Y_test$CAIDE < 4, 0, 1)

roc_list_test <- roc(response = response_test, 
                     predictor = pred_test)

plot(roc_list_test)



