library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(wateRmelon)
library(missMDA)

# Clear workspace and console
rm(list = ls())
cat("\014") 


# Load data
load("~/Data/X_test.RData")
load("~/Data/X_nonTest.RData")
load("~/Data/Y_CAIDE1.RData")
load("~/allModels.RData")

X_all <- cbind(X_test, X_nonTest)
X_all_M <- t(log2(X_all/(1-X_all)))

factors <- names(allModels)
predictedScore_factors <- matrix(NA, nrow = nrow(X_all_M), ncol = length(factors))
colnames(predictedScore_factors) <- factors
rownames(predictedScore_factors) <- rownames(X_all_M)
for (f in 1:length(factors)){
  f1 <- factors[f]
  model <- allModels[[f1]]
  
  # Get features for model fitting
  features <- colnames(model$trainingData)[-10001]
  
  # Make predictions
  predictedScore_factors[,f] <- predict(model, X_all_M[,features], type = "prob")$No
  
}
predictedScore_factors <- as.data.frame(predictedScore_factors)

predictedAge <- agep(X_all, 
                     method='all')

predictedScore_factors$Age <- predictedAge$skinblood.skinblood.age

save(predictedScore_factors, file = "predictedScore_factors_EXTEND.RData")

###############################################################################

# Make predictions

###############################################################################
source("FUN_MachineLearning.R")
load("~/Data/CVindex_CAIDE1.RData")
load("predictedScore_factors_EXTEND.RData")
load("~/Data/Y_CAIDE1.RData")
load("~/Data/Y_test.RData")

X_train <- predictedScore_factors[Y_CAIDE1$Basename,]
Y_train <- Y_CAIDE1
all(rownames(X_train) == Y_train$Basename)


# Settings for repeated cross-valBasename
fitControl <- trainControl(method = "repeatedcv", 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = regressionSummary,
                           index = CVindex)


#*****************************************************************************#
# ElasticNet
#*****************************************************************************#

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


# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$CAIDE,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Save model
save(fit, file = "Fit_CombineFactors_EXTEND_EN.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
X_test <- predictedScore_factors[Y_test$Basename,]
all(rownames(X_test) == Y_test$Basename)
testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE)
RMSE(pred = testPred, obs = Y_test$CAIDE)
R2(pred = testPred, obs = Y_test$CAIDE)

#*****************************************************************************#
# RandomForest
#*****************************************************************************#

library(e1071)
library(ranger)
library(dplyr)

# Number of randomly selected predictors
mtry_CV <- 1:14

# split rule
splitrule_CV <- "variance"

# minimal node size
min.node.size_CV = 1:14

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "ranger"


# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$CAIDE,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Save model
save(fit, file = "Fit_CombineFactors_EXTEND_RF.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
X_test <- predictedScore_factors[Y_test$Basename,]
all(rownames(X_test) == Y_test$Basename)
testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE)
RMSE(pred = testPred, obs = Y_test$CAIDE)
R2(pred = testPred, obs = Y_test$CAIDE)

