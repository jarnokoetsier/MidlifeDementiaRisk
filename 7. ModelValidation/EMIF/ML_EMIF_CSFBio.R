###############################################################################

# Make predictions

###############################################################################

library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library(missMDA)
library(pROC)

# Clear workspace and console
rm(list = ls())
cat("\014") 

load("EMIF/X_train_EMIF.RData")
load("EMIF/Y_train_EMIF.RData")
load("EMIF/X_test_EMIF.RData")
load("EMIF/Y_test_EMIF.RData")

load("EMIF/metaData_fil.RData")
rownames(metaData_fil) <- metaData_fil$X
CSFbio <- metaData_fil[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "Age")]
colnames(CSFbio) <- c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "ChrAge")
samples <- rownames(CSFbio)[(!is.na(CSFbio$Ptau_ASSAY_Zscore)) & 
                              (!is.na(CSFbio$AB_Zscore)) &
                              (!is.na(CSFbio$Ttau_ASSAY_Zscore))]
# MCI and NL only
X_train <- X_train[Y_train$Diagnosis != "AD",]
Y_train <- Y_train[Y_train$Diagnosis != "AD",]
X_test <- X_test[Y_test$Diagnosis != "AD",]
Y_test<- Y_test[Y_test$Diagnosis != "AD",]

#X_train <- X_train[Y_train$Diagnosis == "MCI",]
#Y_train <- Y_train[Y_train$Diagnosis == "MCI",]
#X_test <- X_test[Y_test$Diagnosis == "MCI",]
#Y_test<- Y_test[Y_test$Diagnosis == "MCI",]
#Y_train$Y <- factor(ifelse(Y_train$MCI_Convert == 0 ,"Control","Convert"), levels = c("Control", "Convert"))
#Y_test$Y <- factor(ifelse(Y_test$MCI_Convert == 0 ,"Control","Convert"), levels = c("Control", "Convert"))

Y_train$Y <- factor(ifelse(Y_train$Diagnosis == "NL","Control","MCI"),
                    levels = c("Control", "MCI"))

Y_test$Y <- factor(ifelse(Y_test$Diagnosis == "NL","Control","MCI"),
                   levels = c("Control", "MCI"))


Y_train <- Y_train[intersect(samples, rownames(Y_train)),]
Y_test <- Y_test[intersect(samples, rownames(Y_test)),]

X_train <- cbind.data.frame(X_train[rownames(Y_train),], CSFbio[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], CSFbio[rownames(Y_test),])





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
which(colnames(X_train) == "ChrAge")
set.seed(123)
fit <- train(x = X_train[,-18],
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_CI_EN_CSFbio.RData")

# performance in CV
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha

library(reportROC)
rocdf <- reportROC(gold = Y_test$Y, predictor = testPred$Convert, important = "se")

# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)

# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore", "ChrAge")],
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit_wo, file = "EMIF/Fit_EMIF_CI_EN_CSFbioage.RData")

# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore")],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_wo, file = "EMIF/Fit_EMIF_CI_EN_CSFbioonly.RData")



library(reportROC)
rocdf <- reportROC(gold = Y_test$Y, predictor = testPred$Convert, important = "se")

testPred <- predict(fit_wo, X_test, type = "prob")

roc_test_wo <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_wo)


# Actual training
set.seed(123)
which(colnames(X_train) == "Age")
fit_woa <- train(x = X_train[,-c(14,18)],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_woa, file = "EMIF/Fit_EMIF_CI_EN_CSFbio_noage.RData")

testPred <- predict(fit_woa, X_test, type = "prob")

roc_test_woa <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_woa)



roc.test(roc_test, roc_test_wo)

roc_test <- roc(Y_test$Y, X_test$Age)
auc(roc_test)

roc_test <- roc(Y_test$Y, Y_test$Age)
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
which(colnames(X_train) == "ChrAge")
set.seed(123)
fit <- train(x = X_train[,-18],
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_CI_sPLS_CSFbio.RData")


# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)

# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore", "ChrAge")],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_wo, file = "EMIF/Fit_EMIF_CI_sPLS_CSFbioage.RData")

# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore")],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_wo, file = "EMIF/Fit_EMIF_CI_sPLS_CSFbioonly.RData")

testPred <- predict(fit_wo, X_test, type = "prob")

roc_test_wo <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_wo)

# Actual training
set.seed(123)
fit_woa <- train(x = X_train[,-c(14,18)],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_woa, file = "EMIF/Fit_EMIF_CI_sPLS_CSFbio_noage.RData")

testPred <- predict(fit_woa, X_test, type = "prob")

roc_test_woa <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_woa)



library(reportROC)
rocdf <- reportROC(gold = Y_test$Y, predictor = testPred$MCI, important = "se")



roc.test(roc_test, roc_test_wo)

roc_test <- roc(Y_test$Y, X_test$Age)
auc(roc_test)

roc_test <- roc(Y_test$Y, Y_test$Age)
auc(roc_test)

#*****************************************************************************#
# RandomForest
#*****************************************************************************#
library(e1071)
library(ranger)
library(dplyr)

# Number of randomly selected predictors
mtry_CV <- 1:17

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = 1:17

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "ranger"

# Actual training
which(colnames(X_train) == "ChrAge")
set.seed(123)
fit <- train(x = X_train[,-18],
             y = Y_train$Y,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Save model
save(fit, file = "EMIF/Fit_EMIF_CI_RF_CSFbio.RData")


# Prediction in test set
testPred <- predict(fit, X_test, type = "prob")

roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test)

# Number of randomly selected predictors
mtry_CV <- 1:4

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = 1:4

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "ranger"
# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore", "ChrAge")],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_wo, file = "EMIF/Fit_EMIF_CI_RF_CSFbioage.RData")

# Number of randomly selected predictors
mtry_CV <- 1:3

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = 1:3

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "ranger"
# Actual training
set.seed(123)
fit_wo <- train(x = X_train[,c("Ptau_ASSAY_Zscore", "AB_Zscore", "Ttau_ASSAY_Zscore")],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_wo, file = "EMIF/Fit_EMIF_CI_RF_CSFbioonly.RData")

testPred <- predict(fit_wo, X_test, type = "prob")

roc_test_wo <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_wo)


# Number of randomly selected predictors
mtry_CV <- 1:16

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = 1:16

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "ranger"
set.seed(123)
fit_woa <- train(x = X_train[,-c(14,18)],
                y = Y_train$Y,
                metric= performance_metric,
                method = MLmethod,
                tuneGrid = parameterGrid,
                trControl = fitControl,
                maximize = TRUE)

# Save model
save(fit_woa, file = "EMIF/Fit_EMIF_CI_RF_CSFbio_noage.RData")

testPred <- predict(fit_woa, X_test, type = "prob")

roc_test_woa <- pROC::roc(Y_test$Y, testPred$MCI)
auc(roc_test_woa)



roc.test(roc_test_woa, roc_test_wo)

roc_test <- roc(Y_test$Y, X_test$Age)
auc(roc_test)

roc_test <- roc(Y_test$Y, Y_test$Age)
auc(roc_test)


#*****************************************************************************#
# Plotting
#*****************************************************************************#
load("EMIF/Fit_EMIF_CI_RF_CSFbio.RData")
testPred <- predict(fit, X_test, type = "prob")

load("EMIF/Fit_EMIF_CI_sPLS_CSFbioonly.RData")
testPred_wo <- predict(fit_wo, X_test, type = "prob")

roc_test <- pROC::roc(Y_test$Y, testPred$MCI)
roc_test_wo <- pROC::roc(Y_test$Y, testPred_wo$MCI)
pROC::roc.test(roc_test, roc_test_wo)


load("EMIF/Fit_EMIF_CI_RF_CSFbio.RData")
RF <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_CI_RF_CSFbioonly.RData")
RF_wo <- predict(fit_wo, X_test, type = "prob")
load("EMIF/Fit_EMIF_CI_RF_CSFbio_noage.RData")
RF_woa <- predict(fit_woa, X_test, type = "prob")
load("EMIF/Fit_EMIF_CI_RF_CSFbioage.RData")
RF_wa <- predict(fit_wo, X_test, type = "prob")

load("EMIF/Fit_EMIF_CI_EN_CSFbio.RData")
EN <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_CI_EN_CSFbioonly.RData")
EN_wo <- predict(fit_wo, X_test, type = "prob")
load("EMIF/Fit_EMIF_CI_EN_CSFbio_noage.RData")
EN_woa <- predict(fit_woa, X_test, type = "prob")
load("EMIF/Fit_EMIF_CI_EN_CSFbioage.RData")
EN_wa <- predict(fit_wo, X_test, type = "prob")

load("EMIF/Fit_EMIF_CI_sPLS_CSFbio.RData")
sPLS <- predict(fit, X_test, type = "prob")
load("EMIF/Fit_EMIF_CI_sPLS_CSFbioonly.RData")
sPLS_wo <- predict(fit_wo, X_test, type = "prob")
load("EMIF/Fit_EMIF_CI_sPLS_CSFbio_noage.RData")
sPLS_woa <- predict(fit_woa, X_test, type = "prob")
load("EMIF/Fit_EMIF_CI_sPLS_CSFbioage.RData")
sPLS_wa <- predict(fit_wo, X_test, type = "prob")


testDF <- data.frame(EN = EN$MCI,
                     sPLS = sPLS$MCI,
                     RF = RF$MCI,
                     EN_wo = EN_wo$MCI,
                     sPLS_wo = sPLS_wo$MCI,
                     RF_wo = RF_wo$MCI,
                     EN_woa = EN_woa$MCI,
                     sPLS_woa = sPLS_woa$MCI,
                     RF_woa = RF_woa$MCI,
                     EN_wa = EN_wa$MCI,
                     sPLS_wa = sPLS_wa$MCI,
                     RF_wa = RF_wa$MCI,
                     EpiAge = X_test$Age,
                     Age = Y_test$Age,
                     Sex= Y_test$Gender,
                     Diagnosis = Y_test$Diagnosis,
                     Y = Y_test$Y,
                     Ynum = ifelse(Y_test$Y == "MCI",1,0))



model <- lm(Ynum ~ EN + Age + Sex, data = testDF)
summary(model)

score <- c("EN", "sPLS", "RF")
scoreName <- c("ElasticNet", "sPLS-DA", "Random Forest")
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

plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.17,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- rev(c("#EF3B2C","#CB181D", "#99000D"))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  scale_color_manual(values = colors) +
  ggtitle("CSF biomarkers + MRSs") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_MCI_EMIF_all.png", width = 8, height = 5)




score <- c("EN_wo", "sPLS_wo", "RF_wo")
scoreName <- c("ElasticNet", "sPLS-DA", "Random Forest")
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

plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.17,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- rev(c("#6BAED6","#2171B5","#084594"))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  ggtitle("CSF biomarkers") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_MCI_EMIF_CSFonly.png", width = 8, height = 5)




score <- c("EN_wo", "sPLS_wo", "RF_wo", "EN", "sPLS", "RF")
scoreName <- c("CSF (EN)", "CSF (sPLS-DA)", "CSF (RF)",
               "CSF + MRS (EN)", "CSF + MRS (sPLS-DA)", "CSF + MRS (RF)")
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

plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.34,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- c(rev(c("#6BAED6","#2171B5","#084594")),
            rev(c("#EF3B2C","#CB181D", "#99000D")))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  ggtitle("MCI vs Control") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_MCI_EMIF_combinedPlot.png", width = 8, height = 5)



score <- c("EN_woa", "sPLS_woa", "RF_woa")
scoreName <- c("ElasticNet", "sPLS-DA", "Random Forest")
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

plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.17,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- rev(c("#F16913","#D94801","#8C2D04"))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  ggtitle("CSF biomarkers + MRSs (w/o epi-age)") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_MCI_EMIF_noage.png", width = 8, height = 5)



score <- c("EN_wa", "sPLS_wa", "RF_wa")
scoreName <- c("ElasticNet", "sPLS-DA", "Random Forest")
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

plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.17,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- rev(c("#807DBA","#6A51A3","#4A1486"))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  ggtitle("CSF biomarkers + age") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_MCI_EMIF_age.png", width = 8, height = 5)



load("EMIF/Fit_EMIF_MCI_RF.RData")
testDF$RF1 <- predict(fit, X_test[,1:14], type = "prob")$MCI
load("EMIF/Fit_EMIF_MCI_EN.RData")
testDF$EN1 <- predict(fit, X_test[,1:14], type = "prob")$MCI
load("EMIF/Fit_EMIF_MCI_sPLS.RData")
testDF$sPLS1 <- predict(fit, X_test[,1:14], type = "prob")$MCI


score <- c("EN1", "sPLS1", "RF1", "EN_wo", "sPLS_wo", "RF_wo", "EN", "sPLS", "RF")
scoreName <- c("MRS (EN)", "MRS (sPLS-DA)", "MRS (RF)",
               "CSF (EN)", "CSF (sPLS-DA)", "CSF (RF)",
               "CSF + MRS (EN)", "CSF + MRS (sPLS-DA)", "CSF + MRS (RF)")
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

plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.52,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- c(rev(c("#807DBA","#6A51A3","#4A1486")),
            rev(c("#6BAED6","#2171B5","#084594")),
            rev(c("#EF3B2C","#CB181D", "#99000D")))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_MCI_EMIF_combinedPlot1.png", width = 8, height = 5)




load("EMIF/Fit_EMIF_MCI_RF.RData")
RF1 <- predict(fit, X_test[,1:14], type = "prob")$MCI
test <- pROC::roc(Y_test$Y, RF1)
auc(test)



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

load("EMIF/Fit_EMIF_MCI_RF.RData")
RF1 <- predict(fit, X_test[,1:14], type = "prob")$MCI
test1 <- pROC::roc(Y_test$Y, RF1)
auc(test1)

roc.test(test, test1)
