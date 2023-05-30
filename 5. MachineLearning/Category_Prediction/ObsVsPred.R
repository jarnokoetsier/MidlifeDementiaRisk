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
source("FUN_MachineLearning.R")



###############################################################################

# CAIDE1

###############################################################################

load("Data/Y_CAIDE1.RData")
load("Data/Y_test.RData")
load("Data/CVindex_CAIDE1.RData")
load("CAIDE1_Cor/X_CAIDE1_CorCV.RData")

#****************************************************************************#
# ElasticNet
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"
load(paste0("CV_", Score, "_", FeatureSelection,"_LowRisk_EN.RData"))

X_train = log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))

# Get observed and predicted values
Y_train = Y_CAIDE1
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE < 4, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Set grid for lambda
lambdaCV <- optLambda

# Set grid for alpha
alphaCV <- optAlpha

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]


# Performance on test set
load("CAIDE1_Cor/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE < 4, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$CAIDE

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "CAIDE1_Cor/Performance_CAIDE1_LowRisk_Cor_EN.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"
load(paste0("CV_", Score, "_", FeatureSelection,"_HighRisk_EN.RData"))

X_train = log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))

# Get observed and predicted values
Y_train = Y_CAIDE1
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Set grid for lambda
lambdaCV <- optLambda

# Set grid for alpha
alphaCV <- optAlpha

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]


# Performance on test set
load("CAIDE1_Cor/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE > 7, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$CAIDE

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "CAIDE1_Cor/Performance_CAIDE1_HighRisk_Cor_EN.RData"
)

#****************************************************************************#
# sPLS-DA
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"
load(paste0("CAIDE1_Cor/CV_", Score, "_", FeatureSelection,"_LowRisk_sPLSDA.RData"))

X_train = log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))

# Get observed and predicted values
Y_train = Y_CAIDE1
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE < 4, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Number of component (K)
K_CV <- optK

# Thresholding parameter (eta)
eta_CV <- optEta

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = optKappa

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

MLmethod = "spls"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file ="ObsPred_CV_LowRisk_spls.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]


# Performance on test set
load("CAIDE1_Cor/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE < 4, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$CAIDE

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "CAIDE1_Cor/Performance_CAIDE1_LowRisk_Cor_sPLSDA.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"
load(paste0("CAIDE1_Cor/CV_", Score, "_", FeatureSelection,"_HighRisk_sPLSDA.RData"))

X_train = log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))

# Get observed and predicted values
Y_train = Y_CAIDE1
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Number of component (K)
K_CV <- optK

# Thresholding parameter (eta)
eta_CV <- optEta

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = optKappa

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

MLmethod = "spls"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file ="ObsPred_CV_HighRisk_spls.RData")


# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]


# Performance on test set
load("CAIDE1_Cor/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE > 7, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$CAIDE

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "CAIDE1_Cor/Performance_CAIDE1_HighRisk_Cor_sPLSDA.RData"
)

#****************************************************************************#
# Random forest
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"
load(paste0("CAIDE1_Cor/CV_", Score, "_", FeatureSelection,"_LowRisk_RF.RData"))

X_train = log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))

# Get observed and predicted values
Y_train = Y_CAIDE1
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE < 4, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Number of randomly selected predictors
mtry_CV <- opt_mtry

# split rule
splitrule_CV <- opt_splitrule

# minimal node size
min.node.size_CV = opt_min.node.size

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

MLmethod = "ranger"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file ="ObsPred_CV_LowRisk_RF.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]


# Performance on test set
load("CAIDE1_Cor/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE < 4, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$CAIDE

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "CAIDE1_Cor/Performance_CAIDE1_LowRisk_Cor_RF.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"
load(paste0("CAIDE1_Cor/CV_", Score, "_", FeatureSelection,"_HighRisk_RF.RData"))

X_train = log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))

# Get observed and predicted values
Y_train = Y_CAIDE1
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Number of randomly selected predictors
mtry_CV <- opt_mtry

# split rule
splitrule_CV <- opt_splitrule

# minimal node size
min.node.size_CV = opt_min.node.size

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

MLmethod = "ranger"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file ="ObsPred_CV_HighRisk_RF.RData")


# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]


# Performance on test set
load("CAIDE1_Cor/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE > 7, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$CAIDE

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "CAIDE1_Cor/Performance_CAIDE1_HighRisk_Cor_RF.RData"
)

###############################################################################

# CAIDE2

###############################################################################

load("Data/Y_CAIDE2.RData")
load("Data/Y_test.RData")
load("Data/CVindex_CAIDE2.RData")
load("CAIDE2_Cor/X_CAIDE2_Cor2CV.RData")

#****************************************************************************#
# ElasticNet
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE2"
FeatureSelection = "Cor"
load(paste0("CAIDE2_Cor/CV_", Score, "_", FeatureSelection,"_LowRisk_EN.RData"))

X_train = log2(X_CAIDE2_Cor2CV/(1-X_CAIDE2_Cor2CV))

# Get observed and predicted values
Y_train = Y_CAIDE2
Y_train$Class <- factor(ifelse(Y_CAIDE2$CAIDE < 5, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Set grid for lambda
lambdaCV <- optLambda

# Set grid for alpha
alphaCV <- optAlpha

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
  for (p in 1:8){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1650))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < 8){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file ="ObsPred_CV_LowRisk_en.RData")


# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE2),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_CAIDE2$CAIDE2[unlist(test)]


# Performance on test set
load("CAIDE2_Cor/X_test_Cor2.RData")
testData <- log2(X_test_Cor2/(1-X_test_Cor2))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE2 < 5, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$CAIDE2

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "CAIDE2_Cor/Performance_CAIDE2_LowRisk_Cor_EN.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE2"
FeatureSelection = "Cor"
load(paste0("CAIDE2_Cor/CV_", Score, "_", FeatureSelection,"_HighRisk_EN.RData"))

X_train = log2(X_CAIDE2_Cor2CV/(1-X_CAIDE2_Cor2CV))

# Get observed and predicted values
Y_train = Y_CAIDE2
Y_train$Class <- factor(ifelse(Y_CAIDE2$CAIDE2 > 8, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Set grid for lambda
lambdaCV <- optLambda

# Set grid for alpha
alphaCV <- optAlpha

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
  for (p in 1:8){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1650))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < 8){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file ="ObsPred_CV_HighRisk_en.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE2),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_CAIDE2$CAIDE2[unlist(test)]


# Performance on test set
load("CAIDE2_Cor/X_test_Cor2.RData")
testData <- log2(X_test_Cor2/(1-X_test_Cor2))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE2 > 8, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$CAIDE2

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "CAIDE2_Cor/Performance_CAIDE2_HighRisk_Cor_EN.RData"
)

#****************************************************************************#
# sPLS-DA
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE2"
FeatureSelection = "Cor"
load(paste0("CAIDE2_Cor/CV_", Score, "_", FeatureSelection,"_LowRisk_sPLSDA.RData"))

X_train = log2(X_CAIDE1_Cor2CV/(1-X_CAIDE1_Cor2CV))

# Get observed and predicted values
Y_train = Y_CAIDE2
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE < 5, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Number of component (K)
K_CV <- optK

# Thresholding parameter (eta)
eta_CV <- optEta

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = optKappa

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

MLmethod = "spls"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
  for (p in 1:8){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1650))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < 8){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file ="ObsPred_CV_LowRisk_spls.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE2),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_CAIDE2$CAIDE2[unlist(test)]


# Performance on test set
load("CAIDE1_Cor/X_test_Cor2.RData")
testData <- log2(X_test_Cor2/(1-X_test_Cor2))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE2 < 5, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$CAIDE2

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "CAIDE1_Cor/Performance_CAIDE1_LowRisk_Cor_sPLSDA.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE2"
FeatureSelection = "Cor"
load(paste0("CAIDE2_Cor/CV_", Score, "_", FeatureSelection,"_HighRisk_sPLSDA.RData"))

X_train = log2(X_CAIDE2_Cor2CV/(1-X_CAIDE2_Cor2CV))

# Get observed and predicted values
Y_train = Y_CAIDE2
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 8, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Number of component (K)
K_CV <- optK

# Thresholding parameter (eta)
eta_CV <- optEta

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = optKappa

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

MLmethod = "spls"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
  for (p in 1:8){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1650))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < 8){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file ="ObsPred_CV_HighRisk_spls.RData")


# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE2),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_CAIDE2$CAIDE2[unlist(test)]


# Performance on test set
load("CAIDE2_Cor/X_test_Cor2.RData")
testData <- log2(X_test_Cor2/(1-X_test_Cor2))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE > 8, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$CAIDE2

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "CAIDE2_Cor/Performance_CAIDE2_HighRisk_Cor_sPLSDA.RData"
)

#****************************************************************************#
# Random forest
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE2"
FeatureSelection = "Cor"
load(paste0("CAIDE2_Cor/CV_", Score, "_", FeatureSelection,"_LowRisk_RF.RData"))

X_train = log2(X_CAIDE2_CorCV/(1-X_CAIDE2_CorCV))

# Get observed and predicted values
Y_train = Y_CAIDE1
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE2 < 5, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Number of randomly selected predictors
mtry_CV <- opt_mtry

# split rule
splitrule_CV <- opt_splitrule

# minimal node size
min.node.size_CV = opt_min.node.size

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

MLmethod = "ranger"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
  for (p in 1:8){
    probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1650))
  }
  
  # get exactly 10,000 probes
  n = 1
  finalProbes <- unique(unlist(probes))
  while (length(finalProbes) > 10000){
    probes[[n]] <- probes[[n]][-1]
    finalProbes <- unique(unlist(probes))
    
    if (n < 8){
      n = n + 1
    } else {
      n = 1
    }
  }
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file ="ObsPred_CV_LowRisk_RF.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]


# Performance on test set
load("CAIDE2_Cor/X_test_Cor2.RData")
testData <- log2(X_test_Cor2/(1-X_test_Cor2))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE2 < 5, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$CAIDE2

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "CAIDE1_Cor/Performance_CAIDE2_LowRisk_Cor_RF.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE2"
FeatureSelection = "Cor"
load(paste0("CAIDE2_Cor/CV_", Score, "_", FeatureSelection,"_HighRisk_RF.RData"))

X_train = log2(X_CAIDE2_Cor2CV/(1-X_CAIDE2_Cor2CV))

# Get observed and predicted values
Y_train = Y_CAIDE1
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 8, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Number of randomly selected predictors
mtry_CV <- opt_mtry

# split rule
splitrule_CV <- opt_splitrule

# minimal node size
min.node.size_CV = opt_min.node.size

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

MLmethod = "ranger"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
  for (p in 1:8){
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
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file ="ObsPred_CV_HighRisk_RF.RData")


# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE2),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_CAIDE2$CAIDE2[unlist(test)]


# Performance on test set
load("CAIDE2_Cor/X_test_Cor2.RData")
testData <- log2(X_test_Cor2/(1-X_test_Cor2))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE > 8, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$CAIDE2

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "CAIDE2_Cor/Performance_CAIDE2_HighRisk_Cor_RF.RData"
)

###############################################################################

# LIBRA

###############################################################################

load("Data/Y_LIBRA.RData")
load("Data/Y_test.RData")
load("Data/CVindex_LIBRA.RData")
load("LIBRA_Cor/X_LIBRA_CorLCV.RData")

#****************************************************************************#
# ElasticNet
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "LIBRA"
FeatureSelection = "Cor"
load(paste0(Score, "_Cor/","CV_", Score, "_", FeatureSelection,"_LowRisk_EN.RData"))

X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))

# Get observed and predicted values
Y_train = Y_LIBRA
Y_train$Class <- factor(ifelse(Y_LIBRA$LIBRA < 0, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Set grid for lambda
lambdaCV <- optLambda

# Set grid for alpha
alphaCV <- optAlpha

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file = "ObsPred_CV_lowRisk_LIBRA_EN.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get LIBRA score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_LIBRA),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_LIBRA$LIBRA[unlist(test)]


# Performance on test set
load("LIBRA_Cor/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$LIBRA < 0, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$LIBRA

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "LIBRA_Cor/Performance_LIBRA_LowRisk_Cor_EN.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "LIBRA"
FeatureSelection = "Cor"
load(paste0(Score, "_Cor/","CV_", Score, "_", FeatureSelection,"_HighRisk_EN.RData"))

X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))

# Get observed and predicted values
Y_train = Y_LIBRA
Y_train$Class <- factor(ifelse(Y_LIBRA$LIBRA > 2, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Set grid for lambda
lambdaCV <- optLambda

# Set grid for alpha
alphaCV <- optAlpha

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file = "ObsPred_CV_HighRisk_LIBRA_EN.RData")


# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_LIBRA),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_LIBRA$LIBRA[unlist(test)]


# Performance on test set
load("LIBRA_Cor/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$LIBRA > 2, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$LIBRA

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "LIBRA_Cor/Performance_LIBRA_HighRisk_Cor_EN.RData"
)

#****************************************************************************#
# sPLS-DA
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "LIBRA"
FeatureSelection = "Cor"
load(paste0("LIBRA_Cor/CV_", Score, "_", FeatureSelection,"_LowRisk_sPLSDA.RData"))

X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))

# Get observed and predicted values
Y_train = Y_LIBRA
Y_train$Class <- factor(ifelse(Y_LIBRA$LIBRA < 0, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Number of component (K)
K_CV <- optK

# Thresholding parameter (eta)
eta_CV <- optEta

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = optKappa

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

MLmethod = "spls"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file ="ObsPred_CV_LowRisk_LIBRA_spls.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_LIBRA),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_LIBRA$LIBRA[unlist(test)]


# Performance on test set
load("LIBRA_Cor/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$LIBRA < 0, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$LIBRA

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "LIBRA_Cor/Performance_LIBRA_LowRisk_Cor_sPLSDA.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"
load(paste0("CAIDE1_Cor/CV_", Score, "_", FeatureSelection,"_HighRisk_sPLSDA.RData"))

X_train = log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))

# Get observed and predicted values
Y_train = Y_CAIDE1
Y_train$Class <- factor(ifelse(Y_CAIDE1$CAIDE > 7, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Number of component (K)
K_CV <- optK

# Thresholding parameter (eta)
eta_CV <- optEta

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = optKappa

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

MLmethod = "spls"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file ="ObsPred_CV_HighRisk_spls.RData")


# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]


# Performance on test set
load("CAIDE1_Cor/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$CAIDE > 7, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$CAIDE

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "CAIDE1_Cor/Performance_CAIDE1_HighRisk_Cor_sPLSDA.RData"
)

###############################################################################

# LIBRA

###############################################################################

load("Data/Y_LIBRA.RData")
load("Data/Y_test.RData")
load("Data/CVindex_LIBRA.RData")
load("LIBRA_Cor/X_LIBRA_CorLCV.RData")

#****************************************************************************#
# ElasticNet
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "LIBRA"
FeatureSelection = "Cor"
load(paste0(Score, "_Cor/","CV_", Score, "_", FeatureSelection,"_LowRisk_EN.RData"))

X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))

# Get observed and predicted values
Y_train = Y_LIBRA
Y_train$Class <- factor(ifelse(Y_LIBRA$LIBRA < 0, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Set grid for lambda
lambdaCV <- optLambda

# Set grid for alpha
alphaCV <- optAlpha

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file = "ObsPred_CV_lowRisk_LIBRA_EN.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get LIBRA score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_LIBRA),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_LIBRA$LIBRA[unlist(test)]


# Performance on test set
load("LIBRA_Cor/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$LIBRA < 0, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$LIBRA

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "LIBRA_Cor/Performance_LIBRA_LowRisk_Cor_EN.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "LIBRA"
FeatureSelection = "Cor"
load(paste0(Score, "_Cor/","CV_", Score, "_", FeatureSelection,"_HighRisk_EN.RData"))

X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))

# Get observed and predicted values
Y_train = Y_LIBRA
Y_train$Class <- factor(ifelse(Y_LIBRA$LIBRA > 2, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Set grid for lambda
lambdaCV <- optLambda

# Set grid for alpha
alphaCV <- optAlpha

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

MLmethod = "glmnet"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file = "ObsPred_CV_HighRisk_LIBRA_EN.RData")


# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_LIBRA),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_LIBRA$LIBRA[unlist(test)]


# Performance on test set
load("LIBRA_Cor/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$LIBRA > 2, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$LIBRA

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "LIBRA_Cor/Performance_LIBRA_HighRisk_Cor_EN.RData"
)

#****************************************************************************#
# sPLS-DA
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "LIBRA"
FeatureSelection = "Cor"
load(paste0("LIBRA_Cor/CV_", Score, "_", FeatureSelection,"_LowRisk_sPLSDA.RData"))

X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))

# Get observed and predicted values
Y_train = Y_LIBRA
Y_train$Class <- factor(ifelse(Y_LIBRA$LIBRA < 0, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Number of component (K)
K_CV <- optK

# Thresholding parameter (eta)
eta_CV <- optEta

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = optKappa

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

MLmethod = "spls"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file ="ObsPred_CV_LowRisk_LIBRA_spls.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_LIBRA),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_LIBRA$LIBRA[unlist(test)]


# Performance on test set
load("LIBRA_Cor/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$LIBRA < 0, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$LIBRA

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "LIBRA_Cor/Performance_LIBRA_LowRisk_Cor_sPLSDA.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "LIBRA"
FeatureSelection = "Cor"
load(paste0(Score, "_Cor/","CV_", Score, "_", FeatureSelection,"_HighRisk_sPLSDA.RData"))

X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))

# Get observed and predicted values
Y_train = Y_LIBRA
Y_train$Class <- factor(ifelse(Y_LIBRA$LIBRA > 2, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Number of component (K)
K_CV <- optK

# Thresholding parameter (eta)
eta_CV <- optEta

# kappa (default = 0.5, only relevant for multivariate outcome variables)
kappa_CV = optKappa

# Combine into a single data frame
parameterGrid <- expand.grid(K_CV, eta_CV, kappa_CV)
colnames(parameterGrid) <- c(".K", ".eta", ".kappa")

MLmethod = "spls"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file ="ObsPred_CV_HighRisk_LIBRA_spls.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_LIBRA),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_LIBRA$LIBRA[unlist(test)]


# Performance on test set
load("LIBRA_Cor/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$LIBRA > 2, "High", "Low_Intermediate"),
                                                levels = c("Low_Intermediate", "High")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$LIBRA

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "LIBRA_Cor/Performance_LIBRA_HighRisk_Cor_sPLSDA.RData"
)

#****************************************************************************#
# Random forest
#****************************************************************************#

#============================================================================#
# Low Risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "LIBRA"
FeatureSelection = "Cor"
load(paste0("LIBRA_Cor/CV_", Score, "_", FeatureSelection,"_LowRisk_RF.RData"))

X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))

# Get observed and predicted values
Y_train = Y_LIBRA
Y_train$Class <- factor(ifelse(Y_LIBRA$LIBRA < 0, "Low", "Intermediate_High"),
                        levels = c("Intermediate_High","Low"))

# Number of randomly selected predictors
mtry_CV <- opt_mtry

# split rule
splitrule_CV <- opt_splitrule

# minimal node size
min.node.size_CV = opt_min.node.size

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

MLmethod = "ranger"
performance_metric = "ROC"

ObsPred_CV_LowRisk <- NULL
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
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_LowRisk <- rbind.data.frame(ObsPred_CV_LowRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_LowRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_LowRisk, file ="ObsPred_CV_LowRisk_RF1.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_LowRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_LowRisk$pred <- ifelse(ObsPred_CV_LowRisk$p > threshold, "Intermediate_High","Low")

# Get LIBRRA score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_LIBRA),x)})
ObsPred_CV_LowRisk$ObservedScore <- Y_LIBRA$LIBRA[unlist(test)]


# Performance on test set
load("LIBRA_Cor/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_LowRisk <- data.frame(obs = factor(ifelse(Y_test$LIBRA < 0, "Low", "Intermediate_High"),
                                                levels = c("Intermediate_High", "Low")),
                                   p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_LowRisk$obs, 
                     predictor = ObsPred_test_LowRisk$p)


# Get sensitivies and specificities
plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_LowRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_LowRisk$pred <- ifelse(ObsPred_test_LowRisk$p > threshold, "Intermediate_High", "Low")
ObsPred_test_LowRisk$ObservedScore <- Y_test$LIBRA

# Save
save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "LIBRA_Cor/Performance_LIBRA_LowRisk_Cor_RF.RData"
)

#============================================================================#
# High risk
#============================================================================#

# Load low risk model training
# Score and feature selection method
Score = "LIBRA"
FeatureSelection = "Cor"
load(paste0("LIBRA_Cor/CV_", Score, "_", FeatureSelection,"_HighRisk_RF.RData"))

X_train = log2(X_LIBRA_CorLCV/(1-X_LIBRA_CorLCV))

# Get observed and predicted values
Y_train = Y_LIBRA
Y_train$Class <- factor(ifelse(Y_LIBRA$LIBRA > 2, "High", "Low_Intermediate"),
                        levels = c("Low_Intermediate","High"))

# Number of randomly selected predictors
mtry_CV <- opt_mtry

# split rule
splitrule_CV <- opt_splitrule

# minimal node size
min.node.size_CV = opt_min.node.size

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

MLmethod = "ranger"
performance_metric = "ROC"

ObsPred_CV_HighRisk <- NULL
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
  fitControl <- trainControl(search = "grid", 
                             savePredictions = TRUE,
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
  
  ObsPred_CV_HighRisk <- rbind.data.frame(ObsPred_CV_HighRisk, fit$pred[,1:3])
}
# Set column names
colnames(ObsPred_CV_HighRisk) <- c("pred", "obs", "p")
save(ObsPred_CV_HighRisk, file ="ObsPred_CV_HighRisk_RF.RData")


# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                              Specificity = roc_list_CV$specificities)

# Get AUC
auc_HighRisk_CV <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_HighRisk$pred <- ifelse(ObsPred_CV_HighRisk$p > threshold, "Low_Intermediate","High")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_LIBRA),x)})
ObsPred_CV_HighRisk$ObservedScore <- Y_LIBRA$LIBRA[unlist(test)]


# Performance on test set
load("LIBRA_Cor/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))
pred <- predict(finalModel, t(testData), type = "prob")[,1]
ObsPred_test_HighRisk <- data.frame(obs = factor(ifelse(Y_test$LIBRA > 2, "High", "Low_Intermediate"),
                                                 levels = c("Low_Intermediate", "High")),
                                    p = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test_HighRisk$obs, 
                     predictor = ObsPred_test_HighRisk$p)


# Get sensitivies and specificities
plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                   Specificity = roc_list_test$specificities)

# Get AUC
auc_HighRisk_test <- auc(roc_list_test)

# Get CAIDE1 score
ObsPred_test_HighRisk$pred <- ifelse(ObsPred_test_HighRisk$p > threshold, "Low_Intermediate", "High")
ObsPred_test_HighRisk$ObservedScore <- Y_test$LIBRA

# Save
save(ObsPred_CV_HighRisk,
     auc_HighRisk_CV, 
     threshold,
     ObsPred_test_HighRisk,
     auc_HighRisk_test,
     file = "LIBRA_Cor/Performance_LIBRA_HighRisk_Cor_RF.RData"
)
