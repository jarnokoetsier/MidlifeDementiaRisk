library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(wateRmelon)


###############################################################################

# Make predictions (CAIDE1): PGSs only

###############################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load data
source("FUN_MachineLearning.R")
load("~/PRS/df_list.RData")
load("~/Data/Y_CAIDE1.RData")
load("~/Data/Y_test.RData")

# Prepare data
PGS <- df_list$bayesr
colnames(PGS) <- paste0(colnames(PGS), "_PGS")
samples_train <- Y_CAIDE1[Y_CAIDE1$ID %in% rownames(PGS),]
samples_test <- Y_test[Y_test$ID %in% rownames(PGS),]
X_train <- cbind.data.frame(PGS[samples_train$ID,])
Y_train <- samples_train
all(rownames(X_train) == Y_train$Basename)

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:ncol(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)

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
save(fit, file = "Fit_CombineFactors_CAIDE1_EN_PGSonly.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha

# Prediction in test set
X_test <- cbind.data.frame(PGS[Y_test$ID,])
all(rownames(X_test) == Y_test$ID)

testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE)
RMSE(pred = testPred, obs = Y_test$CAIDE)
R2(pred = testPred, obs = Y_test$CAIDE)

#*****************************************************************************#
# sPLS
#*****************************************************************************#

# Number of component (K)
K_CV <- 1:10

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
save(fit, file = "Fit_CombineFactors_CAIDE1_sPLS_PGSonly.RData")

# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
X_test <- cbind.data.frame(PGS[Y_test$ID,])
all(rownames(X_test) == Y_test$ID)

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
mtry_CV <- seq(6,40,2)

# split rule
splitrule_CV <- "variance"

# minimal node size
min.node.size_CV = 1:15

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
             maximize = FALSE,
             importance = "impurity")

# Save model
save(fit, file = "Fit_CombineFactors_CAIDE1_RF_PGSonly.RData")

# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size



# Prediction in test set
X_test <- cbind.data.frame(PGS[Y_test$ID,])
all(rownames(X_test) == Y_test$ID)

testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE)
RMSE(pred = testPred, obs = Y_test$CAIDE)
R2(pred = testPred, obs = Y_test$CAIDE)


###############################################################################

# Make predictions (CAIDE2): PGSs only

###############################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(wateRmelon)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
source("FUN_MachineLearning.R")
load("~/PRS/df_list.RData")
load("~/Data/Y_CAIDE2.RData")
load("~/Data/Y_test.RData")

# Combine Methylation scores and PGSes
PGS <- df_list$bayesr
colnames(PGS) <- paste0(colnames(PGS), "_PGS")
samples_train <- Y_CAIDE2[Y_CAIDE2$ID %in% rownames(PGS),]
samples_test <- Y_test[Y_test$ID %in% rownames(PGS),]


X_train <- cbind.data.frame(PGS[samples_train$ID,])
Y_train <- samples_train
all(rownames(X_train) == Y_train$Basename)

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:ncol(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)

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
             y = Y_train$CAIDE2,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Save model
save(fit, file = "Fit_CombineFactors_CAIDE2_PGSonly.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
X_test <- cbind.data.frame(PGS[Y_test$ID,])
all(rownames(X_test) == Y_test$ID)

testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE2)
RMSE(pred = testPred, obs = Y_test$CAIDE2)
R2(pred = testPred, obs = Y_test$CAIDE2)

#*****************************************************************************#
# sPLS
#*****************************************************************************#

# Number of component (K)
K_CV <- 1:10

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

# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$CAIDE2,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Save model
save(fit, file = "Fit_CombineFactors_CAIDE2_sPLS_PGSonly.RData")

# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
X_test <- cbind.data.frame(PGS[Y_test$ID,])
all(rownames(X_test) == Y_test$ID)

testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE2)
RMSE(pred = testPred, obs = Y_test$CAIDE2)
R2(pred = testPred, obs = Y_test$CAIDE2)

#*****************************************************************************#
# RandomForest
#*****************************************************************************#

library(e1071)
library(ranger)
library(dplyr)

# Number of randomly selected predictors
mtry_CV <- seq(6,40,2)

# split rule
splitrule_CV <- "variance"

# minimal node size
min.node.size_CV = 1:15

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "ranger"


# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$CAIDE2,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE,
             importance = "impurity")

# Save model
save(fit, file = "Fit_CombineFactors_CAIDE2_RF_PGSonly.RData")

# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size



# Prediction in test set
X_test <- cbind.data.frame(PGS[Y_test$ID,])
all(rownames(X_test) == Y_test$Basename)

testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE2)
RMSE(pred = testPred, obs = Y_test$CAIDE2)
R2(pred = testPred, obs = Y_test$CAIDE2)


###############################################################################

# Make predictions (LIBRA): PGSs only

###############################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(wateRmelon)

# Clear workspace and console
rm(list = ls())
cat("\014") 

source("FUN_MachineLearning.R")
load("predictedScore_factors_EXTEND.RData")
load("~/PRS/df_list.RData")
load("~/Data/Y_LIBRA.RData")
load("~/Data/Y_test.RData")


# Combine Methylation scores and PGSes
PGS <- df_list$bayesr
colnames(PGS) <- paste0(colnames(PGS), "_PGS")
samples_train <- Y_LIBRA[Y_LIBRA$ID %in% rownames(PGS),]
samples_test <- Y_test[Y_test$ID %in% rownames(PGS),]


X_train <- as.data.frame(PGS[samples_train$ID,])
Y_train <- samples_train
all(rownames(X_train) == Y_train$ID)

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:ncol(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)

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
             y = Y_train$LIBRA,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Save model
save(fit, file = "Fit_CombineFactors_LIBRA_EN_PGSonly.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
X_test <- cbind.data.frame(PGS[Y_test$ID,])
all(rownames(X_test) == Y_test$Basename)

testPred <- predict(fit, X_test)

plot(testPred, Y_test$LIBRA)
RMSE(pred = testPred, obs = Y_test$LIBRA)
R2(pred = testPred, obs = Y_test$LIBRA)

#*****************************************************************************#
# sPLS
#*****************************************************************************#

# Number of component (K)
K_CV <- 1:10

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

# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$LIBRA,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Save model
save(fit, file = "Fit_CombineFactors_LIBRA_sPLS_PGSonly.RData")

# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
X_test <- cbind.data.frame(predictedScore_factors[Y_test$Basename,],
                           PGS[Y_test$ID,])
all(rownames(X_test) == Y_test$Basename)

testPred <- predict(fit, X_test)

plot(testPred, Y_test$LIBRA)
RMSE(pred = testPred, obs = Y_test$LIBRA)
R2(pred = testPred, obs = Y_test$LIBRA)

#*****************************************************************************#
# RandomForest
#*****************************************************************************#

library(e1071)
library(ranger)
library(dplyr)

# Number of randomly selected predictors
mtry_CV <- seq(6,40,2)

# split rule
splitrule_CV <- "variance"

# minimal node size
min.node.size_CV = 1:15

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "ranger"


# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$LIBRA,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE,
             importance = "impurity")

# Save model
save(fit, file = "Fit_CombineFactors_LIBRA_RF_PGSonly.RData")

# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size



# Prediction in test set
X_test <- cbind.data.frame(PGS[Y_test$ID,])
all(rownames(X_test) == Y_test$ID)

testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE)
RMSE(pred = testPred, obs = Y_test$LIBRA)
R2(pred = testPred, obs = Y_test$LIBRA)

