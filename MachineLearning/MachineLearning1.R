# Clear workspace and console
rm(list = ls())
cat("\014") 

library(glmnet)
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


#*****************************************************************************#
# Model training CAIDE1
#*****************************************************************************#

# Prepare data (M-values)
X_train = log2(X_CAIDE1_S/(1-X_CAIDE1_S))
Y_train = Y_CAIDE1
all(colnames(X_train) == Y_train$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           number = nfold, 
                           repeats = nrep, 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = regressionSummary)

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0,1,length.out = 11)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "glmnet"

# Register cores for parallel computing
detectCores()
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Actual training
set.seed(123)
fit_CAIDE1 <- train(x = t(X_CAIDE1),
                    y = Y_train$CAIDE,
                    metric= performance_metric,
                    method = MLmethod,
                    tuneGrid = parameterGrid,
                    trControl = fitControl,
                    maximize = FALSE)

# Stop clusters
stopCluster(cl)

# Save model
save(fit_CAIDE1, file = "fit_CAIDE1.RData")
