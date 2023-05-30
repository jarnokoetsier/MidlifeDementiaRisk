################################################################################

# MRS only

################################################################################


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

# Load model
load("EMIF/Fit_EMIF_MCI_sPLS.RData")

# Get predictions
predictions <- fit$pred

# Get threshold
test <- pROC::roc(predictions$obs,predictions$MCI)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))

confusionMatrix(pred,predictions$obs, positive = "MCI")

# Get observed and predicted value
RF <- predict(fit, X_test, type = "prob")
pred <- factor(ifelse(RF$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))
obs <- Y_test$Y
confusionMatrix(obs,pred, positive = "MCI")

library(reportROC)
reportROC(gold = Y_test$Y, predictor = RF$MCI, important = "se")

#==============================================================================#
# AD vs control
#==============================================================================#

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

# MCI and NL only
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

# Load model
load("EMIF/Fit_EMIF_AD_sPLS.RData")

# Get predictions
predictions <- fit$pred

# Get threshold
test <- pROC::roc(predictions$obs,predictions$AD)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$AD > threshold, "AD", "Control"), levels = c("Control", "AD"))

confusionMatrix(pred,predictions$obs, positive = "AD")

# Get observed and predicted value
RF <- predict(fit, X_test, type = "prob")
pred <- factor(ifelse(RF$AD > threshold, "AD", "Control"), levels = c("Control", "AD"))
obs <- Y_test$Y
confusionMatrix(obs,pred, positive = "AD")


library(reportROC)
reportROC(gold = Y_test$Y, predictor = RF$AD, important = "se")

#==============================================================================#
# SCI vs control
#==============================================================================#

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


# Load model
load("EMIF/SCI_prediction/Fit_EMIF_SCIvsControl_sPLS.RData")

# Get predictions
predictions <- fit$pred

# Get threshold
test <- pROC::roc(predictions$obs,predictions$SCI)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$SCI > threshold, "SCI", "Control"), levels = c("NL", "SCI"))

confusionMatrix(pred,predictions$obs, positive = "SCI")

# Get observed and predicted value
RF <- predict(fit, X_test1, type = "prob")
pred <- factor(ifelse(RF$SCI > threshold, "SCI", "NL"), levels = c("NL", "SCI"))
obs <- Y_test1$Y
confusionMatrix(obs,pred, positive = "SCI")

library(reportROC)
reportROC(gold = Y_test1$Y, predictor = 1-RF$SCI, important = "se")

auc(pROC::roc(Y_test1$Y, RF$SCI))

################################################################################

# CSF only

################################################################################

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


# Load model
load("EMIF/Fit_EMIF_CI_RF_CSFbioonly.RData")

# Get predictions
predictions <- fit_wo$pred

# Get threshold
test <- pROC::roc(predictions$obs,predictions$MCI)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))

confusionMatrix(pred,predictions$obs, positive = "MCI")

# Get observed and predicted value
RF <- predict(fit_wo, X_test, type = "prob")
pred <- factor(ifelse(RF$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))
obs <- Y_test$Y
confusionMatrix(obs,pred, positive = "MCI")

library(reportROC)
reportROC(gold = Y_test$Y, predictor = RF$MCI, important = "se")

################################################################################

# CSF + MRS

################################################################################

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


# Load model
load("EMIF/Fit_EMIF_CI_RF_CSFbio.RData")

# Get predictions
predictions <- fit$pred

# Get threshold
test <- pROC::roc(predictions$obs,predictions$MCI)
threshold <- test$thresholds[which.max(sqrt(test$sensitivities*test$specificities))]
pred <- factor(ifelse(predictions$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))

confusionMatrix(pred,predictions$obs, positive = "MCI")

# Get observed and predicted value
RF <- predict(fit, X_test, type = "prob")
pred <- factor(ifelse(RF$MCI > threshold, "MCI", "Control"), levels = c("Control", "MCI"))
obs <- Y_test$Y
confusionMatrix(obs,pred, positive = "MCI")

library(reportROC)
reportROC(gold = Y_test$Y, predictor = RF$MCI, important = "se")
