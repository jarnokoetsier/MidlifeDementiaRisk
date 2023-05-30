# Split data

# Library
library(prospectr)
library(tidyverse)
library(caret)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/metaData_fil.RData")
load("EMIF/predictedScore_factors_fil.RData")

# Male
nTrain0 <- round(0.7*nrow(predictedScore_factors_fil[metaData_fil$Gender == 0,]))
selectedSamples0 <- prospectr::kenStone(
  X = predictedScore_factors_fil[metaData_fil$Gender == 0,], 
  k = nTrain0,
  pc = 8,
  .center = TRUE,
  .scale = TRUE
)
temp0 <- rownames(predictedScore_factors_fil)[metaData_fil$Gender == 0]
test0 <- temp0[selectedSamples0$test]
train0 <- temp0[selectedSamples0$model]

# Female
nTrain1 <- round(0.7*nrow(predictedScore_factors_fil[metaData_fil$Gender == 1,]))
selectedSamples1 <- prospectr::kenStone(
  X = predictedScore_factors_fil[metaData_fil$Gender == 1,], 
  k = nTrain1,
  pc = 8,
  .center = TRUE,
  .scale = TRUE
)
temp1 <- rownames(predictedScore_factors_fil)[metaData_fil$Gender == 1]
test1 <- temp1[selectedSamples1$test]
train1 <- temp1[selectedSamples1$model]


# Combine into single train and test set
X_train <- predictedScore_factors_fil[c(train0, train1),]
Y_train <- metaData_fil[c(train0, train1),]

X_test <- predictedScore_factors_fil[c(test0, test1),]
Y_test <- metaData_fil[c(test0, test1),]

intersect(rownames(X_train), rownames(X_test))

sum(duplicated(X_test))
sum(duplicated(X_train))

table(Y_test$Diagnosis)
table(Y_train$Diagnosis)

save(X_train, file = "EMIF/X_train_EMIF.RData")
save(Y_train, file = "EMIF/Y_train_EMIF.RData")
save(X_test, file = "EMIF/X_test_EMIF.RData")
save(Y_test, file = "EMIF/Y_test_EMIF.RData")

###############################################################################

# Make predictions: NL vs AD

###############################################################################

library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(missMDA)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/X_train_EMIF.RData")
load("EMIF/Y_train_EMIF.RData")
load("EMIF/X_test_EMIF.RData")
load("EMIF/Y_test_EMIF.RData")

# AD and NL only
X_train <- X_train[Y_train$Diagnosis != "MCI",]
Y_train <- Y_train[Y_train$Diagnosis != "MCI",]
X_test <- X_test[Y_test$Diagnosis != "MCI",]
Y_test<- Y_test[Y_test$Diagnosis != "MCI",]

Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","AD"),
                    levels = c("Control", "AD"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","AD"),
                    levels = c("Control", "AD"))

table(Y_test$Y)
table(Y_train$Y)

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)


# Settings for repeated cross-valBasename
fitControl <- trainControl(method = "repeatedcv", 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
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
performance_metric = "ROC"
MLmethod = "glmnet"


# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_AD_EN.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

# AUC
roc_test <- roc(Y_test$Y, testPred$AD)
auc(roc_test)


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
performance_metric = "ROC"
MLmethod = "spls"

# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_AD_sPLS.RData")

# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

# AUC
roc_test <- roc(Y_test$Y, testPred$AD)
auc(roc_test)

#*****************************************************************************#
# RandomForest
#*****************************************************************************#

# Load packages
library(e1071)
library(ranger)
library(dplyr)

# Number of randomly selected predictors
mtry_CV <- 1:14

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = 1:14

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "ranger"

# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE,
             importance = "impurity")

# Save model
save(fit, file = "EMIF/Fit_EMIF_AD_RF.RData")


# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

# AUC
roc_test <- roc(Y_test$Y, testPred$AD)
auc(roc_test)


#*****************************************************************************#
# Evaluate models
#*****************************************************************************#

# Load models
load("EMIF/Fit_EMIF_AD_RF.RData")
RF <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_AD_EN.RData")
EN <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_AD_sPLS.RData")
sPLS <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_CombineFactors_CAIDE1_RF.RData")
EpiCAIDE1 <- predict(fit, X_test)
load("EMIF/Fit_CombineFactors_LIBRA_RF.RData")
EpiLIBRA <- predict(fit, X_test)

# Combine into data frame
testDF <- data.frame(EN = EN$AD,
                     sPLS = sPLS$AD,
                     RF = RF$AD,
                     EpiCAIDE1 = EpiCAIDE1,
                     EpiLIBRA = EpiLIBRA,
                     EpiAge = X_test$Age,
                     Age = Y_test$Age,
                     Sex= Y_test$Gender,
                     Diagnosis = Y_test$Diagnosis,
                     Y = Y_test$Y,
                     Ynum = ifelse(Y_test$Y == "AD",1,0))


# Evaluate significance
model <- lm(Ynum ~ EpiLIBRA + Age + Sex, data = testDF)
summary(model)

# Calculate sensitivities and specificities
score <- c("EN", "sPLS", "RF","EpiCAIDE1", "EpiLIBRA", "EpiAge", "Age")
scoreName <- c("ElasticNet", "sPLS-DA", "Random Forest", "Epi-CAIDE1", "Epi-LIBRA", "Epi-Age", "Age")
plotDF <- as.data.frame(testDF)
ROCplot <- NULL
aucValue <- rep(NA, length(score))
for (i in 1:length(score)){
  test <- pROC::roc(plotDF$Y, plotDF[,score[i]])
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
}

# Format data for plotting
plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.4,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

# Set colors
colors <- rev(c("#E6AB02","#6BAED6","#2171B5","#084594","#EF3B2C","#CB181D", "#99000D"))

# Make ROC curves
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  scale_color_manual(values = colors) +
  ggtitle("AD vs Control") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "EMIF/ROC_AD_EMIF.png", width = 8, height = 5)

###############################################################################

# Make predictions: NL vs MCI

###############################################################################

# Load packages
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
load("EMIF/X_train_EMIF.RData")
load("EMIF/Y_train_EMIF.RData")
load("EMIF/X_test_EMIF.RData")
load("EMIF/Y_test_EMIF.RData")

# MCI and NL only
X_train <- X_train[Y_train$Diagnosis != "AD",]
Y_train <- Y_train[Y_train$Diagnosis != "AD",]
X_test <- X_test[Y_test$Diagnosis != "AD",]
Y_test<- Y_test[Y_test$Diagnosis != "AD",]

Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","MCI"),
                    levels = c("Control", "MCI"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","MCI"),
                   levels = c("Control", "MCI"))

table(Y_test$Y)
table(Y_train$Y)


# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)


# Settings for repeated cross-valBasename
fitControl <- trainControl(method = "repeatedcv", 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
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
performance_metric = "ROC"
MLmethod = "glmnet"


# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCI_EN.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

# AUC
roc_test <- roc(Y_test$Y, testPred$MCI)
auc(roc_test)

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
performance_metric = "ROC"
MLmethod = "spls"

# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCI_sPLS.RData")


# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

# AUC
roc_test <- roc(Y_test$Y, testPred$MCI)
auc(roc_test)


#*****************************************************************************#
# RandomForest
#*****************************************************************************#

# Load packages
library(e1071)
library(ranger)
library(dplyr)

# Number of randomly selected predictors
mtry_CV <- 1:14

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = 1:14

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "ranger"


# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE,
             importance = "impurity")

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCI_RF.RData")


# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size



# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

# AUC
roc_test <- roc(Y_test$Y, testPred$MCI)
auc(roc_test)


#*****************************************************************************#
# Evaluate models
#*****************************************************************************#

# Load models
load("EMIF/Fit_EMIF_MCI_RF.RData")
RF <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_EN.RData")
EN <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_MCI_sPLS.RData")
sPLS <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_CombineFactors_CAIDE1_RF.RData")
EpiCAIDE1 <- predict(fit, X_test)
load("EMIF/Fit_CombineFactors_LIBRA_RF.RData")
EpiLIBRA <- predict(fit, X_test)

# Combine into data frame
testDF <- data.frame(EN = EN$MCI,
                     sPLS = sPLS$MCI,
                     RF = RF$MCI,
                     EpiCAIDE1 = EpiCAIDE1,
                     EpiLIBRA = EpiLIBRA,
                     EpiAge = X_test$Age,
                     Age = Y_test$Age,
                     Sex= Y_test$Gender,
                     Diagnosis = Y_test$Diagnosis,
                     Y = Y_test$Y,
                     Ynum = ifelse(Y_test$Y == "MCI",1,0))

# Evaluate statistical significance
model <- lm(Ynum ~ EpiLIBRA + Age + Sex, data = testDF)
summary(model)

# Calculate sensitivities and specificities
score <- c("EN", "sPLS", "RF","EpiCAIDE1", "EpiLIBRA", "EpiAge", "Age")
scoreName <- c("ElasticNet", "sPLS-DA", "Random Forest", "Epi-CAIDE1", "Epi-LIBRA", "Epi-Age", "Age")
plotDF <- as.data.frame(testDF)
ROCplot <- NULL
aucValue <- rep(NA, length(score))
for (i in 1:length(score)){
  test <- pROC::roc(plotDF$Y, plotDF[,score[i]])
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
}

# Format data for plotting
plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.4,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

# Set colors
colors <- rev(c("#E6AB02","#6BAED6","#2171B5","#084594","#EF3B2C","#CB181D", "#99000D"))

# Make ROC curves
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  scale_color_manual(values = colors) +
  ggtitle("MCI vs Control") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# save plot
ggsave(p, file = "EMIF/ROC_MCI_EMIF1.png", width = 8, height = 5)


