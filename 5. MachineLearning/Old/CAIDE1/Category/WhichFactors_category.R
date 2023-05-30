# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(kernlab)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(pROC)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load cell type composition and phenotype data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Non"

# Load data
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/", f))
}
# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}



methods <- "Non_HighRisk"
load(paste0("CV_CAIDE1/CV_CAIDE1_", methods, ".RData"))

# Performance in test data
testData <- log2(X_test_Non/(1-X_test_Non))
pred_test <- predict(finalModel, t(testData), type = "response")
CAIDE1_cat <- factor(ifelse(Y_test$CAIDE > 7, "High", "Low/Intermediate"),
                          levels = c("Low/Intermediate", "High"))
roc_list_test <- roc(response = CAIDE1_cat, 
                     predictor = pred_test[,1])

perf_test <- as.numeric(auc(roc_list_test))

perf_test_perm <- matrix(NA,100,7)
test_CAIDE1 <- Y_test[,c(1,2,26:32)]
for (factor in 3:9){
  test_CAIDE1_copy <- test_CAIDE1
  set.seed(123)
  for (perm in 1:100){
    rand <- sample(1:nrow(Y_test),nrow(Y_test))
    test_CAIDE1_copy[,factor] <- test_CAIDE1[rand,factor]
    CAIDE1_perm <- rowSums(test_CAIDE1_copy[,3:9])
    CAIDE1_cat_perm <- factor(ifelse(CAIDE1_perm > 7, "High", "Low/Intermediate"),
                              levels = c("Low/Intermediate", "High"))
    
    # Get performance
    roc_list_test <- roc(response = CAIDE1_cat_perm, 
                       predictor = pred_test[,1])
    auc_test <- auc(roc_list_test)
    
    perf_test_perm[perm,factor-2] <- as.numeric(auc_test)
  }
}

colnames(perf_test_perm) <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                              "Total Cholesterol", "Physical Inactivity")
plotDF <- gather(as.data.frame(perf_test_perm))



p <- ggplot(plotDF) +
  geom_hline(yintercept = perf_test, linetype = "dashed",color = "black", linewidth = 1) +
  geom_violin(aes(x = key, y = value, fill = key)) +
  geom_boxplot(aes(x = key, y = value), 
               width=0.1, outlier.shape = NA, position=position_dodge(0.9)) +
  #coord_flip() +
  xlab(NULL) +
  ylab("AUC") +
  labs(fill = NULL) +
  #ggtitle("ElasticNet") +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2") +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5,
                              face = "bold",
                              size = 16))

ggsave(p, file = "WhichFactors_HighRisk_CAIDE1.png", width = 10, height = 6)






methods <- "Non_LowRisk"
load(paste0("CV_CAIDE1/CV_CAIDE1_", methods, ".RData"))

# Performance in test data
testData <- log2(X_test_Non/(1-X_test_Non))
pred_test <- predict(finalModel, t(testData), type = "response")
CAIDE1_cat <- factor(ifelse(Y_test$CAIDE < 4, "Low", "High/Intermediate"),
                     levels = c("High/Intermediate", "Low"))
roc_list_test <- roc(response = CAIDE1_cat, 
                     predictor = pred_test[,1])

perf_test <- as.numeric(auc(roc_list_test))

perf_test_perm <- matrix(NA,100,7)
test_CAIDE1 <- Y_test[,c(1,2,26:32)]
for (factor in 3:9){
  test_CAIDE1_copy <- test_CAIDE1
  set.seed(123)
  for (perm in 1:100){
    rand <- sample(1:nrow(Y_test),nrow(Y_test))
    test_CAIDE1_copy[,factor] <- test_CAIDE1[rand,factor]
    CAIDE1_perm <- rowSums(test_CAIDE1_copy[,3:9])
    CAIDE1_cat_perm <- factor(ifelse(CAIDE1_perm < 4, "Low", "High/Intermediate"),
                              levels = c("High/Intermediate", "Low"))
    
    # Get performance
    roc_list_test <- roc(response = CAIDE1_cat_perm, 
                         predictor = pred_test[,1])
    auc_test <- auc(roc_list_test)
    
    perf_test_perm[perm,factor-2] <- as.numeric(auc_test)
  }
}

colnames(perf_test_perm) <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                              "Total Cholesterol", "Physical Inactivity")
plotDF <- gather(as.data.frame(perf_test_perm))



p <- ggplot(plotDF) +
  geom_hline(yintercept = perf_test, linetype = "dashed",color = "black", linewidth = 1) +
  geom_violin(aes(x = key, y = value, fill = key)) +
  geom_boxplot(aes(x = key, y = value), 
               width=0.1, outlier.shape = NA, position=position_dodge(0.9)) +
  #coord_flip() +
  xlab(NULL) +
  ylab("AUC") +
  labs(fill = NULL) +
  #ggtitle("ElasticNet") +
  theme_minimal() +
  scale_fill_brewer(palette = "Dark2") +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5,
                              face = "bold",
                              size = 16))

ggsave(p, file = "WhichFactors_LowRisk_CAIDE1.png", width = 10, height = 6)