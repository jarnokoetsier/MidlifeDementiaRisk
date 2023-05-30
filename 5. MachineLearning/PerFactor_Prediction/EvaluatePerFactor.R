# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(ranger)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

#*****************************************************************************#
#   AUC of the risk factors
#*****************************************************************************#

# Name of the risk factors
factors <- c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", "Physical", "HeartDisease",
             "Education", "Depression", "SysBP", "Diet")
factorNames <- c("BMI", "Type II Diabetes", "L-M Alchol Intake", "HDL Chol.",
                 "Total Chol.", "Physical Inact.", "Heart Disease", "Education",
                 "Depression", "Systolic BP", "Dietary Intake")

# Get AUC for each model in test set and cross-validation
plotDF_all <- NULL
for (i in 1:length(factors)){
  # Literature-based feature selection + EN
  load("~/PerFactor/Literature/finalOutput_lit_test.RData")
  auc_test_lit <- finalOutput_test[[factors[i]]]$AUC
  
  load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
  auc_CV_lit <- finalOutput[[factors[i]]]$AUC
  
  # Correlation-based feature selection + EN
  load("~/PerFactor/ElasticNet/finalOutput_test.RData")
  auc_test_EN <- finalOutput_test[[factors[i]]]$AUC
  
  load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
  auc_CV_EN <- finalOutput[[factors[i]]]$AUC
  
  # Correlation-based feature selection + RF
  load("~/PerFactor/Random Forest/finalOutput_test.RData")
  auc_test_RF <- finalOutput_test[[factors[i]]]$AUC
  
  load("~/PerFactor/Random Forest/finalOutput_CV.RData")
  auc_CV_RF <- finalOutput[[factors[i]]]$AUC 
  
  
  temp <- data.frame(AUC = c(auc_test_lit, auc_test_EN, auc_test_RF,
                             auc_CV_lit, auc_CV_EN, auc_CV_RF),
                     Factor = rep(factorNames[i],6),
                     Set = c(rep("Test",3),rep("Cross-validation",3)),
                     Method = rep(c("Literature, ElasticNet",
                                    "Correlation, ElasticNet",
                                    "Correlation, Random Forest")))
  
  plotDF_all <- rbind.data.frame(plotDF_all, temp)
}

# Set colors
colors <- c(RColorBrewer::brewer.pal(n = 6,"Dark2"),
            RColorBrewer::brewer.pal(n = 6,"Set2"))

# Make plots
p <- ggplot(plotDF_all) +
  geom_bar(aes(x = Factor, y = AUC, fill = Factor, alpha = Method),
           stat = "identity", position = position_dodge()) +
  facet_grid(rows = vars(Set)) +
  coord_cartesian(ylim = c(0.5,1)) +
  scale_alpha_manual(values = c(0.5,0.75,1)) +
  scale_fill_manual(values = colors) +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")

# Save plots
ggsave(p, file = "AUC_all_Methylation.png", width = 8, height = 6)


#*****************************************************************************#
#   Case-control distribution
#*****************************************************************************#

# Test set distribution
n_cases <- rep(NA, length(finalOutput_test))
n_control <- rep(NA, length(finalOutput_test))
for (i in 1:length(finalOutput_test)){
  n_cases[i] <- sum(finalOutput_test[[i]]$ObsPred_test$obs == "Yes")
  n_control[i] <- sum(finalOutput_test[[i]]$ObsPred_test$obs == "No")
}

plotDF_n <- data.frame(Number = c(n_cases,n_control),
                       CaseControl = c(rep("Case", length(n_cases)),
                                       rep("Control", length(n_control))),
                       Factor = factorNames)

# Make plot
p <- ggplot(plotDF_n)+
  geom_bar(aes(x = Factor, y = Number, fill = CaseControl), 
           stat = "identity", color = "black") +
  coord_flip() +
  ylab("Number of samples") +
  xlab(NULL) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

# Save plot
ggsave(p, file = "PerFactor/CaseControlDistribution.png", width = 8, height = 6)

# Cross-validatio distribution
load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
n_cases <- rep(NA, length(finalOutput))
n_control <- rep(NA, length(finalOutput))
for (i in 1:length(finalOutput)){
  n_cases[i] <- sum(finalOutput[[i]]$ObsPred_CV$obs == "Yes")/5
  n_control[i] <- sum(finalOutput[[i]]$ObsPred_CV$obs == "No")/5
}

plotDF_n <- data.frame(Number = c(n_cases,n_control,924-n_cases-n_control),
                       CaseControl = c(rep("Case", length(n_cases)),
                                       rep("Control", length(n_control)),
                                       rep("NA", length(n_control))),
                       Factor = factorNames)

# Set colors
colors <- c(RColorBrewer::brewer.pal(n = 3, "Set1")[1:2], "#595959")

# Make plot
p <- ggplot(plotDF_n)+
  geom_bar(aes(x = Factor, y = Number, fill = CaseControl), 
           stat = "identity", color = "black") +
  coord_flip() +
  ylab("Number of samples") +
  xlab(NULL) +
  theme_minimal() +
  scale_fill_manual(values = colors) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

# Save plot
ggsave(p, file = "PerFactor/CaseControlDistribution_training.png", width = 8, height = 6)