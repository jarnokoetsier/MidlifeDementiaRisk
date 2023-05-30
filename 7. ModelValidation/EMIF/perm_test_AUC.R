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

X_train <- X_train[Y_train$Diagnosis != "AD",]
Y_train <- Y_train[Y_train$Diagnosis != "AD",]
X_test <- X_test[Y_test$Diagnosis != "AD",]
Y_test<- Y_test[Y_test$Diagnosis != "AD",]

load("EMIF/Fit_EMIF_MCI_RF.RData")
RF1 <- predict(fit, X_test, type = "prob")$MCI
test <- pROC::roc(factor(Y_test$Diagnosis, levels = c("NL", "MCI")), RF1)
original_perf <- as.numeric(auc(test))



load("EMIF/metaData_fil.RData")
rownames(metaData_fil) <- metaData_fil$X
CSFbio <- metaData_fil[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "Age")]
colnames(CSFbio) <- c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "ChrAge")
samples <- rownames(CSFbio)[(!is.na(CSFbio$Ptau_ASSAY_Zscore)) & 
                              (!is.na(CSFbio$AB_Zscore)) &
                              (!is.na(CSFbio$Ttau_ASSAY_Zscore))]
samples <- intersect(samples, rownames(Y_test))


RF1 <- predict(fit, X_test[samples,], type = "prob")$MCI
test <- pROC::roc(factor(Y_test[samples,"Diagnosis"], levels = c("NL", "MCI")), RF1)
selected_perf <- as.numeric(auc(test))



table(Y_test$Diagnosis)
table(Y_test[samples,"Diagnosis"])
n <- 10000
perm_perf <- rep(NA,n)
set.seed(123)
for (i in 1:n){
  # Randomly select 85 MCI
  perm_samples_MCI <- rownames(Y_test)[Y_test$Diagnosis == "MCI"][sample(1:88,85)]
  
  # Randomly select 20 NL
  perm_samples_NL <- rownames(Y_test)[Y_test$Diagnosis == "NL"][sample(1:72,20)]
  
  perm_samples <- c(perm_samples_MCI,perm_samples_NL)
  
  RF1 <- predict(fit, X_test[perm_samples,], type = "prob")$MCI
  test <- pROC::roc(factor(Y_test[perm_samples,"Diagnosis"], levels = c("NL", "MCI")), RF1, quiet = TRUE)
  perm_perf[i] <- as.numeric(auc(test))
  
}

sum(perm_perf > selected_perf)

plotDF <- data.frame(Value = perm_perf,
                     ID = 1:length(perm_perf))

p <- ggplot(plotDF) +
  geom_histogram(aes(x = Value), fill = "grey") +
  geom_vline(xintercept = original_perf) +
  geom_label(aes(x = original_perf, y = 100, label = "Complete test set\nperformance")) +
  geom_vline(xintercept = selected_perf) +
  geom_label(aes(x = selected_perf, y = 100, label = "CSF test set\nperformance")) +
  xlab("AUC") +
  ylab("Count") +
  theme_bw()

ggsave(p, file = "AUC_permutation.png", width = 8, height = 5)