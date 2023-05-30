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
methods = c("cor", "Cor_spls", "Cor_svm")
methodNames <- c("ElasticNet", "sPLS", "SVM")
perf_test <- rep(NA, length(methods))
plotDF_all <- NULL
for (m in 1:length(methods)){
  load(paste0("CV_CAIDE1/CV_CAIDE1_", methods[m], ".RData"))
  
  # Performance in test data
  testData <- log2(X_test_Cor/(1-X_test_Cor))
  pred_test <- predict(finalModel, t(testData))
  perf_test[m] <- RMSE(pred = pred_test,
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
  plotDF$Method <- rep(methodNames[m],nrow(plotDF))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

p <- ggplot(plotDF_all) +
  geom_hline(yintercept = perf_test[1], linetype = "dashed",color = "#FD841F", linewidth = 1) +
  geom_hline(yintercept = perf_test[2], linetype = "dashed",color = "#E14D2A", linewidth = 1) +
  geom_hline(yintercept = perf_test[3], linetype = "dashed",color = "#CD104D", linewidth = 1) +
  geom_violin(aes(x = key, y = value, fill = Method)) +
  geom_boxplot(aes(x = key, y = value, group = interaction(key,Method)), 
               width=0.1, outlier.shape = NA, position=position_dodge(0.9)) +
  #coord_flip() +
  xlab(NULL) +
  ylab("RMSE") +
  labs(fill = NULL) +
  #ggtitle("ElasticNet") +
  theme_minimal() +
  scale_fill_manual(values = c("#FD841F", "#E14D2A", "#CD104D")) +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

ggsave(p, file = "ImportanceFactors_CAIDE1_all_horizontal.png", width = 10, height = 6)



################################################################################

# Cohen's f2

################################################################################

#Score and feature selection method
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
methods = c("cor", "Cor_spls", "Cor_svm")
methodNames <- c("ElasticNet", "sPLS", "SVM")
plotDF_all <- NULL
for (m in 1:length(methods)){
  load(paste0("CV_CAIDE1/CV_CAIDE1_", methods[m], ".RData"))
  
  # Performance in test data
  testData <- log2(X_test_Cor/(1-X_test_Cor))
  pred_test <- predict(finalModel, t(testData))
  
  dataMatrix <- cbind.data.frame(Y_test[,c(26:32)],pred_test)
  colnames(dataMatrix) <- c(colnames(Y_test[,c(26:32)]), "PredictedScore")
  
  formula <- paste0("PredictedScore ~ ", 
                    paste0("0 + ", paste(colnames(Y_test[,26:32]), collapse = " + ")))
  
                            
  model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
  
  # Get fitted and residual values
  fittedValues <- fitted(model)
  residualValues <- residuals(model)
  
  # Calculate R-squared
  sse <- sum(residualValues^2)
  ssr <- sum(fittedValues^2)
  Rsquared = ssr/(ssr + sse)
  
  # Get global effect size
  globalEffect <- (Rsquared)/(1-Rsquared)
  
  cohenF <- list()
  Y_factors <- Y_test[,26:32]
  for (factor in 1:7){
    formula <- paste0("PredictedScore ~ ", 
                      paste0("0 + ", paste(colnames(Y_factors[,-factor]), collapse = " + ")))
    
    
    model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
    
    # Get fitted and residual values
    fittedValues <- fitted(model)
    residualValues <- residuals(model)
    
    # Calculate R-squared
    sse <- sum(residualValues^2)
    ssr <- sum(fittedValues^2)
    Rsquared_factor = ssr/(ssr + sse)
    
    # Calculate cohen's f2 statistic (local effect size)
    cohenF[[factor]] <- ((Rsquared - Rsquared_factor)/(1-Rsquared))
  }
  
  factorNames <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                                "Total Cholesterol", "Physical Inactivity")
  plotDF <- data.frame(Name = c(factorNames, "Global"), 
                       cohenF2  = c(unlist(cohenF), globalEffect))
  plotDF$Method <- rep(methodNames[m],nrow(plotDF))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}


p <- ggplot(plotDF_all[plotDF$Name != "Global",])+
  geom_bar(aes(x = Name, y = cohenF2, fill = Method),
           stat = "identity", position=position_dodge()) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression("Cohen's " ~ f^2)) +
  labs(fill = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#FD841F", "#E14D2A", "#CD104D")) +
  theme(legend.position = "bottom")

ggsave(p, file = "WhichFactors_CohenF2.png", width = 8, height = 6)

















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
