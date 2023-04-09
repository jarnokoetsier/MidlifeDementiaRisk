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

# Make predictions (CAIDE1)

###############################################################################

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
save(fit, file = "Fit_CombineFactors_CAIDE1_sPLS.RData")

# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

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
             maximize = FALSE,
             importance = "impurity")

# Save model
save(fit, file = "Fit_CombineFactors_EXTEND_RF.RData")

# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size



# Prediction in test set
X_test <- predictedScore_factors[Y_test$Basename,]
all(rownames(X_test) == Y_test$Basename)
testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE)
RMSE(pred = testPred, obs = Y_test$CAIDE)
R2(pred = testPred, obs = Y_test$CAIDE)
cor.test(testPred, Y_test$CAIDE, method = "spearman")


###############################################################################

# Make predictions (CAIDE2)

###############################################################################

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

source("FUN_MachineLearning.R")
load("~/Data/CVindex_CAIDE2.RData")
load("predictedScore_factors_EXTEND.RData")
load("~/Data/Y_CAIDE2.RData")
load("~/Data/Y_test.RData")

X_train <- predictedScore_factors[Y_CAIDE2$Basename,]
Y_train <- Y_CAIDE2
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
             y = Y_train$CAIDE2,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Save model
save(fit, file = "Fit_CombineFactors_CAIDE2_EN.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
X_test <- predictedScore_factors[Y_test$Basename,]
all(rownames(X_test) == Y_test$Basename)
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
save(fit, file = "Fit_CombineFactors_CAIDE2_sPLS.RData")

# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
X_test <- predictedScore_factors[Y_test$Basename,]
all(rownames(X_test) == Y_test$Basename)
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
             y = Y_train$CAIDE2,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE,
             importance = "impurity")

# Save model
save(fit, file = "Fit_CombineFactors_CAIDE2_RF.RData")

# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size



# Prediction in test set
X_test <- predictedScore_factors[Y_test$Basename,]
all(rownames(X_test) == Y_test$Basename)
testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE2)
RMSE(pred = testPred, obs = Y_test$CAIDE2)
R2(pred = testPred, obs = Y_test$CAIDE2)
cor.test(testPred, Y_test$CAIDE2, method = "spearman")


###############################################################################

# Make predictions (LIBRA)

###############################################################################
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

source("FUN_MachineLearning.R")
load("~/Data/CVindex_LIBRA.RData")
load("predictedScore_factors_EXTEND.RData")
load("~/Data/Y_LIBRA.RData")
load("~/Data/Y_test.RData")

X_train <- predictedScore_factors[Y_LIBRA$Basename,]
Y_train <- Y_LIBRA
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
             y = Y_train$LIBRA,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Save model
save(fit, file = "Fit_CombineFactors_LIBRA_EN.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
X_test <- predictedScore_factors[Y_test$Basename,]
all(rownames(X_test) == Y_test$Basename)
testPred <- predict(fit, X_test)

plot(testPred, Y_test$LIBRA)
RMSE(pred = testPred, obs = Y_test$LIBRA)
R2(pred = testPred, obs = Y_test$LIBRA)
cor.test(testPred, Y_test$LIBRA, method = "spearman")

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
save(fit, file = "Fit_CombineFactors_LIBRA_sPLS.RData")

# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
X_test <- predictedScore_factors[Y_test$Basename,]
all(rownames(X_test) == Y_test$Basename)
testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE)
RMSE(pred = testPred, obs = Y_test$LIBRA)
R2(pred = testPred, obs = Y_test$LIBRA)
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
             y = Y_train$LIBRA,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE,
             importance = "impurity")

# Save model
save(fit, file = "Fit_CombineFactors_LIBRA_RF.RData")

# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size



# Prediction in test set
X_test <- predictedScore_factors[Y_test$Basename,]
all(rownames(X_test) == Y_test$Basename)
testPred <- predict(fit, X_test)

plot(testPred, Y_test$CAIDE)
RMSE(pred = testPred, obs = Y_test$LIBRA)
R2(pred = testPred, obs = Y_test$LIBRA)


factorNames <- c("BMI", "Type II Diabetes", "L-M Alcohol Intake", "HDL Cholesterol",
                 "Total Cholesterol", "Physical Inactivity", "Heart Disease", "Education",
                 "Depression", "Systolic Blood Pressure", "Dietary Intake", "Sex", "Smoking",
                 "Age")

# Variable importance
load("~/Fit_CombineFactors_LIBRA_RF.RData")
plotDF_LIBRA <- data.frame(Gini = varImp(fit, scale = TRUE)$importance,
                           Factor = factorNames,
                           Method = rep("LIBRA", length(varImp(fit)$importance)))

load("~/Fit_CombineFactors_CAIDE1_RF.RData")
plotDF_CAIDE1 <- data.frame(Gini = varImp(fit, scale = TRUE)$importance,
                           Factor = factorNames,
                           Method = rep("CAIDE1", length(varImp(fit)$importance)))

load("~/Fit_CombineFactors_CAIDE2_RF.RData")
plotDF_CAIDE2 <- data.frame(Gini = varImp(fit, scale = TRUE)$importance,
                            Factor = factorNames,
                            Method = rep("CAIDE2", length(varImp(fit)$importance)))

plotDF_all <- rbind.data.frame(plotDF_CAIDE1, plotDF_CAIDE2, plotDF_LIBRA)

p <- ggplot(plotDF_all) +
  geom_bar(aes(x = Factor, y = Overall, fill = Method),
           stat = "identity",
           position = position_dodge()) +
  xlab(NULL) +
  ylab("Scaled Gini Index") +
  scale_fill_manual(values = c("#FB6A4A","#FEC44F", "#BCBDDC")) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

ggsave(p, file = "varImp_RF_perFactor.png", width = 8, height = 6)
