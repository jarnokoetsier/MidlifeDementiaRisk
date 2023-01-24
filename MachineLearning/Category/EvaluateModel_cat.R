rm(list = ls())
cat("\014") 
gc()

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Load machine learning functions
source("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA/MachineLearning/FUN_MachineLearning.R")

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Low Risk
load("CV_CAIDE1/Performance_CAIDE1_LowRisk_Cor_EN.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_CV_LowRisk$p)

# Get sensitivies and specificities
plotDF_LowRisk_CV <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)


roc_list_test <- roc(response = factor(ObsPred_test_LowRisk$obs,
                                     levels = c("Intermediate_High", "Low")), 
                   predictor = ObsPred_test_LowRisk$p)

plotDF_LowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                             Specificity = roc_list_test$specificities)

# High Risk
load("CV_CAIDE1/Performance_CAIDE1_HighRisk_Cor_EN.RData")

# construct ROC
roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                     levels = c("Low_Intermediate", "High")), 
                   predictor = ObsPred_CV_HighRisk$p)

# Get sensitivies and specificities
plotDF_HighRisk_CV <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                                Specificity = roc_list_CV$specificities)


roc_list_test <- roc(response = factor(ObsPred_test_HighRisk$obs,
                                       levels = c("Low_Intermediate", "High")), 
                     predictor = ObsPred_test_HighRisk$p)

plotDF_HighRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

library(RColorBrewer)
colors <- c(brewer.pal(n = 4, name = "Reds")[c(2,4)],
            brewer.pal(n = 4, name = "Blues")[c(2,4)])

p <- ggplot() +
  geom_path(data = plotDF_LowRisk_CV, aes(y = Sensitivity, x = 1- Specificity,
                                          color = "Low Risk (CV)"), 
            linewidth = 2, linetype = "solid") +
  geom_path(data = plotDF_LowRisk_test, aes(y = Sensitivity, x = 1- Specificity,
                                            color = "Low Risk (Test)"), 
            linewidth = 2, linetype = "solid") +
  geom_path(data = plotDF_HighRisk_CV, aes(y = Sensitivity, x = 1- Specificity,
                                         color = "High Risk (CV)"), 
           linewidth = 2, linetype = "solid") +
  geom_path(data = plotDF_HighRisk_test, aes(y = Sensitivity, x = 1- Specificity,
                                            color = "High Risk (Test)"),
            linewidth = 2, linetype = "solid") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  scale_color_manual(values = colors) +
  ggtitle("CAIDE1",
          subtitle = "ElasticNet-regularized logistic regression model") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))


ggsave(p, file = "CAIDE1_Cat_Cor_AUC_EN.png", height = 8, width = 10)




load("CV_CAIDE1/Performance_CAIDE1_HighRisk_Cor_EN.RData")
# G

ObsPred_CV_LowRisk$pred <- str_replace_all(ObsPred_CV_LowRisk$pred,"/", "_")
ObsPred_CV_LowRisk$obs <- str_replace_all(ObsPred_CV_LowRisk$obs,"/", "_")
ObsPred_test_LowRisk$pred <- str_replace_all(ObsPred_test_LowRisk$pred,"/", "_")
ObsPred_test_LowRisk$obs <- str_replace_all(ObsPred_test_LowRisk$obs,"/", "_")

save(ObsPred_CV_LowRisk,
     auc_LowRisk_CV, 
     threshold,
     ObsPred_test_LowRisk,
     auc_LowRisk_test,
     file = "Performance_CAIDE1_LowRisk_Cor_EN.RData"
)

###############################################################################

# Test set

###############################################################################

# Load test data
load("X/X_Cor_CAIDE1/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))

method <- c("EN")
methodName <- c("ElasticNet")

auc_lowRisk_test <- NULL
plotDF_lowRisk_test <- NULL
auc_highRisk_test <- NULL
plotDF_highRisk_test <- NULL
for (i in 1:length(method)){
  
  #****************************************************************************#
  # Low Risk Model
  #****************************************************************************#
  
  # Load data
  load(paste0("CV_CAIDE1/CV_CAIDE1_Cor_LowRisk_", method[i],".RData"))
  
  # Get predicted value
  pred <- predict(finalModel, t(testData), type = "prob")
  
  # Combine predicted and observed into data frame
  ObsPred_test <- data.frame(Observed = ifelse(Y_test$CAIDE < 4, "Low", "Intermediate/High"),
                             Predicted = as.numeric(pred$Low))
  
  # Calculate AUC
  roc_list_test <- roc(response = factor(ObsPred_test$Observed, levels = c("Intermediate/High", "Low")), 
                       predictor = ObsPred_test$Predicted)
  
  auc_lowRisk_test_temp <- data.frame(AUC = as.numeric(auc(roc_list_test)),
                                 Method = methodName[i])
  
  auc_lowRisk_test <- rbind.data.frame(auc_lowRisk_test, auc_lowRisk_test_temp)
  
  
  # Get sensitivies and specificities
  plotDF_lowRisk_test_temp <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                    Specificity = roc_list_test$specificities,
                                    Method = rep(methodName[i], length(roc_list_test$specificities)))
  
  plotDF_lowRisk_test <- rbind.data.frame(plotDF_lowRisk_test, plotDF_lowRisk_test_temp)
  
  
  #****************************************************************************#
  # High Risk Model
  #****************************************************************************#
  
  # Load data
  load(paste0("CV_CAIDE1/CV_CAIDE1_Cor_highRisk_", method[i],".RData"))
  
  # Get predicted value
  pred <- predict(finalModel, t(testData), type = "prob")
  
  # Combine predicted and observed into data frame
  ObsPred_test <- data.frame(Observed = ifelse(Y_test$CAIDE > 7, "High", "Low/Intermediate"),
                             Predicted = as.numeric(pred$High))
  
  # Calculate AUC
  roc_list_test <- roc(response = factor(ObsPred_test$Observed, levels = c("Low/Intermediate", "High")), 
                       predictor = ObsPred_test$Predicted)
  
  auc_highRisk_test_temp <- data.frame(AUC = as.numeric(auc(roc_list_test)),
                                      Method = methodName[i])
  
  auc_highRisk_test <- rbind.data.frame(auc_highRisk_test, auc_highRisk_test_temp)
  
  
  # Get sensitivies and specificities
  plotDF_highRisk_test_temp <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                         Specificity = roc_list_test$specificities,
                                         Method = rep(methodName[i], length(roc_list_test$specificities)))
  
  plotDF_highRisk_test <- rbind.data.frame(plotDF_highRisk_test, plotDF_highRisk_test_temp)
  
  
}

p <- ggplot() +
  geom_path(data = plotDF_lowRisk_test, aes(y = Sensitivity, x = 1- Specificity,
                                       color = Method), 
            linewidth = 2) +
  geom_path(data = plotDF_highRisk_test, aes(y = Sensitivity, x = 1- Specificity,
                                        color = "High Risk (8-14)"), 
            linewidth = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  ggtitle("CAIDE1") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) +
  scale_color_brewer(palette = "Set1")

ggsave(p, file = "CAIDE1_Cat_Cor_AUC_lowRisk.png", height = 8, width = 8)



p <- ggplot() +

  geom_path(data = plotDF_highRisk_test, aes(y = Sensitivity, x = 1- Specificity,
                                             color = Method), 
            linewidth = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  ggtitle("CAIDE1") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) +
  scale_color_brewer(palette = "Set1")

ggsave(p, file = "CAIDE1_Cat_Cor_AUC_highRisk.png", height = 8, width = 8)

###############################################################################

# cross-validation

###############################################################################

