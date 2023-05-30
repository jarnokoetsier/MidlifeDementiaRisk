library(tidyverse)
library(caret)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/metaData_EMIF.RData")
load("EMIF/predictedScore_factors_EMIF.RData")

metaData_EMIF <- as.data.frame(metaData_EMIF)
rownames(metaData_EMIF) <- metaData_EMIF$X

# Remove individuals with unmatching sex
metaData_EMIF <- metaData_EMIF[metaData_EMIF$sex.match == 1,]

# Keep midlife samples only
metaData_EMIF <- metaData_EMIF[metaData_EMIF$Age <= 75,]

# Remove converters
converters1 <-  unique(rownames(metaData_EMIF)[metaData_EMIF$CTR_Convert == 1])[-1]
converters1 <- intersect(converters1, metaData_EMIF$X[(metaData_EMIF$Diagnosis == "NL")])
converters2 <- unique(metaData_EMIF$X[(metaData_EMIF$LastFU_Diagnosis == "MCI") | (metaData_EMIF$LastFU_Diagnosis == "AD")])[-1]
converters2 <- intersect(converters2, metaData_EMIF$X[(metaData_EMIF$Diagnosis == "NL")])
converters_all <- unique(c(converters1, converters2))
metaData_EMIF <- metaData_EMIF[setdiff(rownames(metaData_EMIF), converters_all),]
table(metaData_EMIF$CTR_Convert)
table(metaData_EMIF$Diagnosis)

# Get SCI
SCI <- unique(rownames(metaData_EMIF)[metaData_EMIF$Diagnosis == "SCI"])
X_SCI <- predictedScore_factors[SCI,]
Y_SCI <- metaData_EMIF[SCI,]

# Male
nTrain0 <- 38
selectedSamples0 <- prospectr::kenStone(
  X = X_SCI, 
  k = nTrain0,
  pc = 3,
  .center = TRUE,
  .scale = TRUE
)

# Combine into single train and test set
X_train_SCI <- X_SCI[selectedSamples0$model,]
Y_train_SCI <- Y_SCI[selectedSamples0$model,]

X_test_SCI <- X_SCI[selectedSamples0$test,]
Y_test_SCI <- Y_SCI[selectedSamples0$test,]

intersect(rownames(X_train_SCI), rownames(X_test_SCI))

sum(duplicated(X_test_SCI))
sum(duplicated(X_train_SCI))

table(Y_test_SCI$Diagnosis)
table(Y_train_SCI$Diagnosis)

# Save SCI data
save(X_train_SCI, file = "EMIF/X_train_SCI_EMIF.RData")
save(Y_train_SCI, file = "EMIF/Y_train_SCI_EMIF.RData")
save(X_test_SCI, file = "EMIF/X_test_SCI_EMIF.RData")
save(Y_test_SCI, file = "EMIF/Y_test_SCI_EMIF.RData")

###############################################################################

# Make predictions: MCI vs SCI

###############################################################################

# Load packages
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
load("EMIF/X_train_SCI_EMIF.RData")
load("EMIF/Y_train_SCI_EMIF.RData")
load("EMIF/X_test_SCI_EMIF.RData")
load("EMIF/Y_test_SCI_EMIF.RData")

# AD and NL only
X_train1 <- rbind.data.frame(X_train[Y_train$Diagnosis == "MCI",], X_train_SCI)
Y_train1 <-  rbind.data.frame(Y_train[Y_train$Diagnosis == "MCI",], Y_train_SCI)
X_test1 <-  rbind.data.frame(X_test[Y_test$Diagnosis == "MCI",], X_test_SCI)
Y_test1 <- rbind.data.frame(Y_test[Y_test$Diagnosis == "MCI",], Y_test_SCI)

Y_train1$Y <- factor(Y_train1$Diagnosis,
                    levels = c("SCI", "MCI"))

Y_test1$Y <- factor(Y_test1$Diagnosis,
                   levels = c("SCI", "MCI"))

table(Y_test1$Y)
table(Y_train1$Y)

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train1),5, returnTrain = TRUE)
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
fit <- train(x = X_train1,
             y = Y_train1$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCIvsSCI_EN.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
testPred <- predict(fit, X_test1, type = "prob")
roc_test <- pROC::roc(Y_test1$Y, testPred$MCI)
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
fit <- train(x = X_train1,
             y = Y_train1$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCIvsSCI_sPLS.RData")


# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
testPred <- predict(fit, X_test1, type = "prob")
roc_test <- pROC::roc(Y_test1$Y, testPred$MCI)
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
fit <- train(x = X_train1,
             y = Y_train1$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE,
             importance = "impurity")

# Save model
save(fit, file = "EMIF/Fit_EMIF_MCIvsSCI_RF.RData")


# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size

# Prediction in test set
testPred <- predict(fit, X_test1, type = "prob")

# AUC
roc_test <- pROC::roc(Y_test1$Y, testPred$MCI)
auc(roc_test)


#*****************************************************************************#
# Evaluate models
#*****************************************************************************#

# Load models
load("EMIF/SCI_prediction/Fit_EMIF_MCIvsSCI_RF.RData")
RF <- predict(fit, X_test1, type = "prob")
load("EMIF/SCI_prediction/Fit_EMIF_MCIvsSCI_EN.RData")
EN <- predict(fit, X_test1, type = "prob")
load("EMIF/SCI_prediction/Fit_EMIF_MCIvsSCI_sPLS.RData")
sPLS <- predict(fit, X_test1, type = "prob")
load("EMIF/Fit_CombineFactors_CAIDE1_RF.RData")
EpiCAIDE1 <- predict(fit, X_test1)
load("EMIF/Fit_CombineFactors_LIBRA_RF.RData")
EpiLIBRA <- predict(fit, X_test1)

# Combine into data frame
testDF <- data.frame(EN = EN$MCI,
                     sPLS = sPLS$MCI,
                     RF = RF$MCI,
                     EpiCAIDE1 = EpiCAIDE1,
                     EpiLIBRA = EpiLIBRA,
                     EpiAge = X_test1$Age,
                     Age = Y_test1$Age,
                     Sex= Y_test1$Gender,
                     Diagnosis = Y_test1$Diagnosis,
                     Y = Y_test1$Y,
                     Ynum = ifelse(Y_test1$Y == "MCI",1,0))

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
  ggtitle("MCI vs SCI") +
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
ggsave(p, file = "EMIF/ROC_MCIvsSCI_test_EMIF.png", width = 7.5, height = 5)

###############################################################################

# Make predictions: SCI vs Normal

###############################################################################

# Load packages
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
load("EMIF/X_train_SCI_EMIF.RData")
load("EMIF/Y_train_SCI_EMIF.RData")
load("EMIF/X_test_SCI_EMIF.RData")
load("EMIF/Y_test_SCI_EMIF.RData")

# AD and NL only
X_train1 <- rbind.data.frame(X_train[Y_train$Diagnosis == "NL",], X_train_SCI)
Y_train1 <-  rbind.data.frame(Y_train[Y_train$Diagnosis == "NL",], Y_train_SCI)
X_test1 <-  rbind.data.frame(X_test[Y_test$Diagnosis == "NL",], X_test_SCI)
Y_test1 <- rbind.data.frame(Y_test[Y_test$Diagnosis == "NL",], Y_test_SCI)

Y_train1$Y <- factor(Y_train1$Diagnosis,
                     levels = c("NL", "SCI"))

Y_test1$Y <- factor(Y_test1$Diagnosis,
                    levels = c("NL", "SCI"))

table(Y_test1$Y)
table(Y_train1$Y)

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train1),5, returnTrain = TRUE)
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
fit <- train(x = X_train1,
             y = Y_train1$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_SCIvsControl_EN.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
testPred <- predict(fit, X_test1, type = "prob")

roc_test <- pROC::roc(Y_test1$Y, testPred$SCI)
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
fit <- train(x = X_train1,
             y = Y_train1$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_SCIvsControl_sPLS.RData")


# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
testPred <- predict(fit, X_test1, type = "prob")

roc_test <- pROC::roc(Y_test1$Y, testPred$SCI)
auc(roc_test)


#*****************************************************************************#
# RandomForest
#*****************************************************************************#

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
fit <- train(x = X_train1,
             y = Y_train1$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE,
             importance = "impurity")

# Save model
save(fit, file = "EMIF/Fit_EMIF_SCIvsControl_RF.RData")


# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size



# Prediction in test set
testPred <- predict(fit, X_test1, type = "prob")

# AUC
roc_test <- pROC::roc(Y_test1$Y, testPred$SCI)
auc(roc_test)


#*****************************************************************************#
# Evaluate models
#*****************************************************************************#

# Load models
load("EMIF/SCI_prediction/Fit_EMIF_SCIvsControl_RF.RData")
RF <- predict(fit, X_test1, type = "prob")
load("EMIF/SCI_prediction/Fit_EMIF_SCIvsControl_EN.RData")
EN <- predict(fit, X_test1, type = "prob")
load("EMIF/SCI_prediction/Fit_EMIF_SCIvsControl_sPLS.RData")
sPLS <- predict(fit, X_test1, type = "prob")
load("EMIF/Fit_CombineFactors_CAIDE1_RF.RData")
EpiCAIDE1 <- predict(fit, X_test1)
load("EMIF/Fit_CombineFactors_LIBRA_RF.RData")
EpiLIBRA <- predict(fit, X_test1)

# Combine into data frame
testDF <- data.frame(EN = EN$SCI,
                     sPLS = sPLS$SCI,
                     RF = RF$SCI,
                     EpiCAIDE1 = EpiCAIDE1,
                     EpiLIBRA = EpiLIBRA,
                     EpiAge = X_test1$Age,
                     Age = Y_test1$Age,
                     Sex= Y_test1$Gender,
                     Diagnosis = Y_test1$Diagnosis,
                     Y = Y_test1$Y,
                     Ynum = ifelse(Y_test1$Y == "SCI",1,0))

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
  ggtitle("SCI vs Control") +
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
ggsave(p, file = "EMIF/ROC_SCIvsControl_test_EMIF.png", width = 7.5, height = 5)


###############################################################################

# Make predictions: AD vs SCI

###############################################################################

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(missMDA)

# Clear workspace and console
rm(list = ls())
cat("\014") 

load("EMIF/X_train_EMIF.RData")
load("EMIF/Y_train_EMIF.RData")
load("EMIF/X_test_EMIF.RData")
load("EMIF/Y_test_EMIF.RData")
load("EMIF/X_train_SCI_EMIF.RData")
load("EMIF/Y_train_SCI_EMIF.RData")
load("EMIF/X_test_SCI_EMIF.RData")
load("EMIF/Y_test_SCI_EMIF.RData")

# AD and NL only
X_train1 <- rbind.data.frame(X_train[Y_train$Diagnosis == "AD",], X_train_SCI)
Y_train1 <-  rbind.data.frame(Y_train[Y_train$Diagnosis == "AD",], Y_train_SCI)
X_test1 <-  rbind.data.frame(X_test[Y_test$Diagnosis == "AD",], X_test_SCI)
Y_test1 <- rbind.data.frame(Y_test[Y_test$Diagnosis == "AD",], Y_test_SCI)

Y_train1$Y <- factor(Y_train1$Diagnosis,
                     levels = c("SCI", "AD"))

Y_test1$Y <- factor(Y_test1$Diagnosis,
                    levels = c("SCI", "AD"))

table(Y_test1$Y)
table(Y_train1$Y)

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train1),5, returnTrain = TRUE)
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
fit <- train(x = X_train1,
             y = Y_train1$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/SCI_prediction/Fit_EMIF_ADvsSCI_EN.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha


# Prediction in test set
testPred <- predict(fit, X_test1, type = "prob")

roc_test <- pROC::roc(Y_test1$Y, testPred$AD)
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
fit <- train(x = X_train1,
             y = Y_train1$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/SCI_prediction/Fit_EMIF_ADvsSCI_sPLS.RData")


# performance in CV
trainResults <- fit$results
optK <- fit$bestTune$K
optEta <- fit$bestTune$eta
optKappa <- fit$bestTune$kappa

# Prediction in test set
testPred <- predict(fit, X_test1, type = "prob")

roc_test <- pROC::roc(Y_test1$Y, testPred$AD)
auc(roc_test)


#*****************************************************************************#
# RandomForest
#*****************************************************************************#

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
fit <- train(x = X_train1,
             y = Y_train1$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE,
             importance = "impurity")

# Save model
save(fit, file = "EMIF/SCI_prediction/Fit_EMIF_ADvsSCI_RF.RData")


# performance in CV
trainResults <- fit$results
opt_mtry <- fit$bestTune$mtry
opt_splitrule <- fit$bestTune$splitrule
opt_min.node.size = fit$bestTune$min.node.size

# Prediction in test set
testPred <- predict(fit, X_test1, type = "prob")

# AUC
roc_test <- pROC::roc(Y_test1$Y, testPred$AD)
auc(roc_test)

#*****************************************************************************#
# Evaluate models
#*****************************************************************************#

# Load models
load("EMIF/SCI_prediction/Fit_EMIF_ADvsSCI_RF.RData")
RF <- predict(fit, X_test1, type = "prob")
load("EMIF/SCI_prediction/Fit_EMIF_ADvsSCI_EN.RData")
EN <- predict(fit, X_test1, type = "prob")
load("EMIF/SCI_prediction/Fit_EMIF_ADvsSCI_sPLS.RData")
sPLS <- predict(fit, X_test1, type = "prob")
load("EMIF/Fit_CombineFactors_CAIDE1_RF.RData")
EpiCAIDE1 <- predict(fit, X_test1)
load("EMIF/Fit_CombineFactors_LIBRA_RF.RData")
EpiLIBRA <- predict(fit, X_test1)

# Combine into data frame
testDF <- data.frame(EN = EN$AD,
                     sPLS = sPLS$AD,
                     RF = RF$AD,
                     EpiCAIDE1 = EpiCAIDE1,
                     EpiLIBRA = EpiLIBRA,
                     EpiAge = X_test1$Age,
                     Age = Y_test1$Age,
                     Sex= Y_test1$Gender,
                     Diagnosis = Y_test1$Diagnosis,
                     Y = Y_test1$Y,
                     Ynum = ifelse(Y_test1$Y == "AD",1,0))

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
  ggtitle("MCI vs SCI") +
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
ggsave(p, file = "EMIF/ROC_ADvsSCI_test_EMIF.png", width = 7.5, height = 5)




