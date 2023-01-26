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
library(RColorBrewer)

# Load machine learning functions
source("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA/MachineLearning/FUN_MachineLearning.R")

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}


################################################################################

# ROC

################################################################################

methods <- c("EN", "sPLSDA")
methodNames <- c("ElasticNet",  "sPLS-DA")

# Low Risk
plotDF_LowRisk_CV <- NULL
plotDF_LowRisk_test <- NULL
for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_CAIDE1/Performance_CAIDE1_LowRisk_Cor_", methods[i],".RData"))
  
  # construct ROC
  roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                       levels = c("Intermediate_High", "Low")), 
                     predictor = ObsPred_CV_LowRisk$p)
  
  # Get sensitivies and specificities
  temp <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                     Specificity = roc_list_CV$specificities,
                     Method = paste0(rep(methodNames[i], length(roc_list_CV$specificities))," (Low Risk)"))
  
  plotDF_LowRisk_CV <- rbind.data.frame(plotDF_LowRisk_CV, temp)
  
  
  roc_list_test <- roc(response = factor(ObsPred_test_LowRisk$obs,
                                         levels = c("Intermediate_High", "Low")), 
                       predictor = ObsPred_test_LowRisk$p)
  
  temp <- data.frame(Sensitivity = roc_list_test$sensitivities,
                     Specificity = roc_list_test$specificities,
                     Method = paste0(rep(methodNames[i], length(roc_list_test$specificities))," (Low Risk)"))
  
  plotDF_LowRisk_test <- rbind.data.frame(plotDF_LowRisk_test, temp)
  
}

# High Risk
plotDF_HighRisk_CV <- NULL
plotDF_HighRisk_test <- NULL
for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_CAIDE1/Performance_CAIDE1_HighRisk_Cor_", methods[i],".RData"))
  
  # construct ROC
  roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                       levels = c("Low_Intermediate", "High")), 
                     predictor = ObsPred_CV_HighRisk$p)
  
  # Get sensitivies and specificities
  temp <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                     Specificity = roc_list_CV$specificities,
                     Method = paste0(rep(methodNames[i], length(roc_list_CV$specificities))," (High Risk)"))
  
  plotDF_HighRisk_CV <- rbind.data.frame(plotDF_HighRisk_CV, temp)
  
  
  roc_list_test <- roc(response = factor(ObsPred_test_HighRisk$obs,
                                         levels = c("Low_Intermediate", "High")), 
                       predictor = ObsPred_test_HighRisk$p)
  
  temp <- data.frame(Sensitivity = roc_list_test$sensitivities,
                     Specificity = roc_list_test$specificities,
                     Method = paste0(rep(methodNames[i], length(roc_list_test$specificities))," (High Risk)"))
  
  plotDF_HighRisk_test <- rbind.data.frame(plotDF_HighRisk_test, temp)
  
}



colors <- c(brewer.pal(n = 4, name = "Reds")[c(3)],
            brewer.pal(n = 4, name = "Blues")[c(3)],
            brewer.pal(n = 4, name = "Reds")[c(4)],
            brewer.pal(n = 4, name = "Blues")[c(4)])

p <- ggplot() +
  geom_path(data = plotDF_LowRisk_CV, aes(y = Sensitivity, x = 1- Specificity,
                                          color = Method), 
            linewidth = 1.5, linetype = "solid") +
  geom_path(data = plotDF_HighRisk_CV, aes(y = Sensitivity, x = 1- Specificity,
                                            color = Method), 
            linewidth = 1.5, linetype = "solid") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  scale_color_manual(values = colors) +
  ggtitle("CAIDE1",
          subtitle = "Performance in cross-validation") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))


ggsave(p, file = "CAIDE1_Cat_Cor_AUC_CV.png", height = 6, width = 9)


p <- ggplot() +
  geom_path(data = plotDF_LowRisk_test, aes(y = Sensitivity, x = 1- Specificity,
                                          color = Method), 
            linewidth = 1.5, linetype = "solid") +
  geom_path(data = plotDF_HighRisk_test, aes(y = Sensitivity, x = 1- Specificity,
                                           color = Method), 
            linewidth = 1.5, linetype = "solid") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  scale_color_manual(values = colors) +
  ggtitle("CAIDE1",
          subtitle = "Performance in test set") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))


ggsave(p, file = "CAIDE1_Cat_Cor_AUC_test.png", height = 6, width = 9)


################################################################################

# Boxplots

################################################################################
methods <- c("EN", "sPLSDA")
methodNames <- c("ElasticNet",  "sPLS-DA")

ObsPred_CV_all <- NULL
ObsPred_test_all <- NULL
for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_CAIDE1/Performance_CAIDE1_HighRisk_Cor_", methods[i],".RData"))
  load(paste0("CV_CAIDE1/Performance_CAIDE1_LowRisk_Cor_", methods[i],".RData"))
  
  # Obs vs predicted in CV
  temp <- data.frame(Score = ObsPred_CV_HighRisk$ObservedScore,
                               Class = rep("Intermediate Risk (4-7)",nrow(ObsPred_CV_HighRisk)))
  temp$Class[(ObsPred_CV_HighRisk$pred == "High") &
                         (ObsPred_CV_LowRisk$pred == "Intermediate_High")] <- "High Risk (8-14)"
  temp$Class[(ObsPred_CV_LowRisk$pred == "Low") &
                         (ObsPred_CV_HighRisk$pred == "Low_Intermediate")] <- "Low Risk (0-3)"
  temp$Method <- rep(methodNames[i], nrow(temp))
  
  ObsPred_CV_all <- rbind.data.frame(ObsPred_CV_all, temp)
  
  # Obs vs predicted in test set
  temp <- data.frame(Score = ObsPred_test_HighRisk$ObservedScore,
                     Class = rep("Intermediate Risk (4-7)",nrow(ObsPred_test_HighRisk)))
  temp$Class[(ObsPred_test_HighRisk$pred == "High") &
               (ObsPred_test_LowRisk$pred == "Intermediate_High")] <- "High Risk (8-14)"
  temp$Class[(ObsPred_test_LowRisk$pred == "Low") &
               (ObsPred_test_HighRisk$pred == "Low_Intermediate")] <- "Low Risk (0-3)"
  temp$Method <- rep(methodNames[i], nrow(temp))
  
  ObsPred_test_all <- rbind.data.frame(ObsPred_test_all, temp)
  
}


ObsPred_CV_all$Class <- factor(ObsPred_CV_all$Class,
                               levels = c("Low Risk (0-3)","Intermediate Risk (4-7)","High Risk (8-14)"))
p <- ggplot(ObsPred_CV_all) +
  geom_rect(ymin = -Inf, ymax = 4, xmin = -Inf, xmax = Inf,fill = "#FFF5F0", alpha = 0.5) +
  geom_rect(ymin = 4, ymax = 7, xmin = -Inf, xmax = Inf,fill = "#FEE0D2", alpha = 0.5) +
  geom_rect(ymin = 7, ymax = Inf, xmin = -Inf, xmax = Inf,fill = "#FCBBA1", alpha = 0.5) +
  geom_hline(yintercept = 4, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 7, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_boxplot(aes(x = Class, y = Score, fill = Method), 
               alpha = 1,outlier.shape = NA) +
  geom_point(aes(x= Class, y = Score, color = Method, group = Method),
             position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0.2)) +
  xlab("Predicted Class") +
  ylab("CAIDE1 Score") +
  ggtitle("CAIDE1",subtitle = "Performance in cross-validation") +
  scale_fill_manual(values = c("#FB6A4A","#6BAED6")) +
  scale_color_manual(values = c("#A50F15","#08519C")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(p, file = "CAIDE1_Cat_Cor_boxPlot_CV.png", height = 6, width = 9)



ObsPred_test_all$Class <- factor(ObsPred_test_all$Class,
                               levels = c("Low Risk (0-3)","Intermediate Risk (4-7)","High Risk (8-14)"))
p <- ggplot(ObsPred_test_all) +
  geom_rect(ymin = -Inf, ymax = 4, xmin = -Inf, xmax = Inf,fill = "#FFF5F0", alpha = 0.5) +
  geom_rect(ymin = 4, ymax = 7, xmin = -Inf, xmax = Inf,fill = "#FEE0D2", alpha = 0.5) +
  geom_rect(ymin = 7, ymax = Inf, xmin = -Inf, xmax = Inf,fill = "#FCBBA1", alpha = 0.5) +
  geom_hline(yintercept = 4, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 7, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_boxplot(aes(x = Class, y = Score, fill = Method), 
               alpha = 1,outlier.shape = NA) +
  geom_point(aes(x= Class, y = Score, color = Method, group = Method),
             position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0.2)) +
  xlab("Predicted Class") +
  ylab("CAIDE1 Score") +
  ggtitle("CAIDE1",subtitle = "Performance in test set") +
  scale_fill_manual(values = c("#FB6A4A","#6BAED6")) +
  scale_color_manual(values = c("#A50F15","#08519C")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(p, file = "CAIDE1_Cat_Cor_boxPlot_CV.png", height = 6, width = 9)


