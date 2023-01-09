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

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load cell type composition and phenotype data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

################################################################################

# ElasticNet

################################################################################

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"

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

# Load data
methods = "cor"
load(paste0("CV_CAIDE1/CV_CAIDE1_", methods, ".RData"))

# Performance in test data
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred_test <- predict(finalModel, t(testData))
perf_test <- RMSE(pred = pred_test,
                  obs = Y_test$CAIDE)

perf_test_perm <- matrix(NA,100,7)
test_CAIDE1 <- Y_test[,c(1,2,26:32)]
for (factor in 3:9){
  test_CAIDE1_copy <- test_CAIDE1
  set.seed(123)
  for (perm in 1:100){
    rand <- sample(1:nrow(Y_test),nrow(Y_test))
    test_CAIDE1_copy[,factor] <- test_CAIDE1[rand,factor]
    CAIDE1_perm <- rowSums(test_CAIDE1_copy[,3:9])
    
    perf_test_perm[perm,factor-2] <- RMSE(pred = pred_test,
                                          obs = CAIDE1_perm)
  }
}

colnames(perf_test_perm) <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                              "Total Cholesterol", "Physical Inactivity")
plotDF <- gather(as.data.frame(perf_test_perm))


p <- ggplot(plotDF) +
  geom_hline(yintercept = perf_test, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_violin(aes(x = key, y = value, fill = key)) +
  geom_boxplot(aes(x = key, y = value), width = 0.1) +
  xlab("CAIDE1 Factors") +
  ylab("RMSE") +
  ggtitle("ElasticNet") +
  theme_classic() +
  #scale_fill_hue(c = 100, l = 40) +
  scale_fill_brewer(palette = "Accent") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

ggsave(p, file = "ImportanceFactors_CAIDE1_EN.png", width = 8, height = 6)

################################################################################

# sPLS

################################################################################

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"

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

# Load data
methods = "Cor_spls"
load(paste0("CV_CAIDE1/CV_CAIDE1_", methods, ".RData"))

# Performance in test data
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred_test <- predict(finalModel, t(testData))
perf_test <- RMSE(pred = pred_test,
                  obs = Y_test$CAIDE)

perf_test_perm <- matrix(NA,100,7)
test_CAIDE1 <- Y_test[,c(1,2,26:32)]
for (factor in 3:9){
  test_CAIDE1_copy <- test_CAIDE1
  set.seed(123)
  for (perm in 1:100){
    rand <- sample(1:nrow(Y_test),nrow(Y_test))
    test_CAIDE1_copy[,factor] <- test_CAIDE1[rand,factor]
    CAIDE1_perm <- rowSums(test_CAIDE1_copy[,3:9])
    
    perf_test_perm[perm,factor-2] <- RMSE(pred = pred_test,
                                          obs = CAIDE1_perm)
  }
}

colnames(perf_test_perm) <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                              "Total Cholesterol", "Physical Inactivity")
plotDF <- gather(as.data.frame(perf_test_perm))


p <- ggplot(plotDF) +
  geom_hline(yintercept = perf_test, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_violin(aes(x = key, y = value, fill = key)) +
  geom_boxplot(aes(x = key, y = value), width = 0.1) +
  xlab("CAIDE1 Factors") +
  ylab("RMSE") +
  ggtitle("sPLS") +
  theme_classic() +
  #scale_fill_hue(c = 100, l = 40) +
  scale_fill_brewer(palette = "Accent") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

ggsave(p, file = "ImportanceFactors_CAIDE1_sPLS.png", width = 8, height = 6)


################################################################################

# SVM

################################################################################

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"

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

# Load data
methods = "Cor_svm"
load(paste0("CV_CAIDE1/CV_CAIDE1_", methods, ".RData"))

# Performance in test data
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred_test <- predict(finalModel, t(testData))
perf_test <- RMSE(pred = pred_test,
                  obs = Y_test$CAIDE)

perf_test_perm <- matrix(NA,100,7)
test_CAIDE1 <- Y_test[,c(1,2,26:32)]
for (factor in 3:9){
  test_CAIDE1_copy <- test_CAIDE1
  set.seed(123)
  for (perm in 1:100){
    rand <- sample(1:nrow(Y_test),nrow(Y_test))
    test_CAIDE1_copy[,factor] <- test_CAIDE1[rand,factor]
    CAIDE1_perm <- rowSums(test_CAIDE1_copy[,3:9])
    
    perf_test_perm[perm,factor-2] <- RMSE(pred = pred_test,
                                          obs = CAIDE1_perm)
  }
}

colnames(perf_test_perm) <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                              "Total Cholesterol", "Physical Inactivity")
plotDF <- gather(as.data.frame(perf_test_perm))


p <- ggplot(plotDF) +
  geom_hline(yintercept = perf_test, linetype = "dashed", color = "black", linewidth = 1.5) +
  geom_violin(aes(x = key, y = value, fill = key)) +
  geom_boxplot(aes(x = key, y = value), width = 0.1) +
  xlab("CAIDE1 Factors") +
  ylab("RMSE") +
  ggtitle("SVM") +
  theme_classic() +
  #scale_fill_hue(c = 100, l = 40) +
  scale_fill_brewer(palette = "Accent") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

ggsave(p, file = "ImportanceFactors_CAIDE1_SVM.png", width = 8, height = 6)
