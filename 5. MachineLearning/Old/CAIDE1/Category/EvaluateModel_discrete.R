# Clear workspace and console
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

# Load data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Non"

# Load CV index
load("CVindex_CAIDE1.RData")


################################################################################

# ROC in CV

################################################################################

load("X/X_Non/X_test_Non.RData")
testData <-  log2(X_test_Non/(1-X_test_Non))

#****************************************************************************#
# Cross-validation and test set
#****************************************************************************#

#============================================================================#
# Low risk
#============================================================================#
# Load low risk model training
load(paste0("CV_CAIDE1/CV_", Score, "_", FeatureSelection,"_LowRisk.RData"))

# construct ROC
roc_list_CV <- roc(response = ObsPred_CV$Observed, 
                   predictor = ObsPred_CV$Predicted)

# Get sensitivies and specificities
plotDF_lowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)

# Get AUC
auc_lowRisk <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_lowRisk <- ObsPred_CV
ObsPred_CV_lowRisk$PredictedClass <- ifelse(ObsPred_CV_lowRisk$Predicted > threshold, "Low", "High/Intermediate")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_lowRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]


# Performance on test set
pred <- predict(finalModel, t(testData), type = "response")
ObsPred_test <- data.frame(Observed = ifelse(Y_test$CAIDE < 4, 2, 1),
                           Predicted = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test$Observed, 
                   predictor = ObsPred_test$Predicted)


# Get sensitivies and specificities
plotDF_lowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_lowRisk_test <- auc(roc_list_test)

ObsPred_test_lowRisk <- ObsPred_test
ObsPred_test_lowRisk$PredictedClass <- ifelse(ObsPred_test_lowRisk$Predicted > threshold, "Low", "High/Intermediate")
ObsPred_test_lowRisk$ObservedScore <- Y_test$CAIDE

#============================================================================#
# High risk
#============================================================================#

# Load high risk model
load(paste0("CV_CAIDE1/CV_", Score, "_", FeatureSelection,"_HighRisk.RData"))

# Construct ROC
roc_list_CV <- roc(response = ObsPred_CV$Observed, 
                   predictor = ObsPred_CV$Predicted)

# Get sensitivies and specificities
plotDF_highRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                             Specificity = roc_list_CV$specificities)
# Get AUC
auc_highRisk <- auc(roc_list_CV)

# Get predictions
Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
ObsPred_CV_highRisk <- ObsPred_CV
ObsPred_CV_highRisk$PredictedClass <- ifelse(ObsPred_CV_highRisk$Predicted > threshold, "High", "Low/Intermediate")

# Get CAIDE score
test <- lapply(CVindex,function(x){setdiff(1:nrow(Y_CAIDE1),x)})
ObsPred_CV_highRisk$ObservedScore <- Y_CAIDE1$CAIDE[unlist(test)]

# Performance on test set
pred <- predict(finalModel, t(testData), type = "response")
ObsPred_test <- data.frame(Observed = ifelse(Y_test$CAIDE > 7, 2, 1),
                           Predicted = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test$Observed, 
                     predictor = ObsPred_test$Predicted)


# Get sensitivies and specificities
plotDF_highRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_highRisk_test <- auc(roc_list_test)

ObsPred_test_highRisk <- ObsPred_test
ObsPred_test_highRisk$PredictedClass <- ifelse(ObsPred_test_highRisk$Predicted > threshold, "High", "Low/Intermediate")
ObsPred_test_highRisk$ObservedScore <- Y_test$CAIDE

#============================================================================#
# plot ROC
#============================================================================#

p <- ggplot() +
  geom_path(data = plotDF_lowRisk, aes(y = Sensitivity, x = 1- Specificity,
                                       color = "Low Risk (1-3)"), 
            linewidth = 2) +
  geom_path(data = plotDF_highRisk, aes(y = Sensitivity, x = 1- Specificity,
                                        color = "High Risk (8-14)"), 
            linewidth = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_text(aes(x = 0.75, y = 0.25, label = paste0("AUC (Low Risk): ", round(auc_lowRisk,2), "\n"),
                color = "Low Risk (1-3)"), fontface = "bold") +
  geom_text(aes(x = 0.75, y = 0.25, label = paste0("\nAUC (High Risk): ", round(auc_highRisk,2)),
                color = "High Risk (8-14)"), fontface = "bold") +
  ggtitle("Performance in Cross-Validation") +
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

ggsave(p, file = "CAIDE1_Non_EN_Discrete_AUC.png", height = 8, width = 8)


p <- ggplot() +
  geom_path(data = plotDF_lowRisk_test, aes(y = Sensitivity, x = 1- Specificity,
                                       color = "Low Risk (1-3)"), 
            linewidth = 2) +
  geom_path(data = plotDF_highRisk_test, aes(y = Sensitivity, x = 1- Specificity,
                                        color = "High Risk (8-14)"), 
            linewidth = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_text(aes(x = 0.75, y = 0.25, label = paste0("AUC (Low Risk): ", round(auc_lowRisk_test,2), "\n"),
                color = "Low Risk (1-3)"), fontface = "bold") +
  geom_text(aes(x = 0.75, y = 0.25, label = paste0("\nAUC (High Risk): ", round(auc_highRisk_test,2)),
                color = "High Risk (8-14)"), fontface = "bold") +
  ggtitle("Performance in Test Set") +
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

ggsave(p, file = "CAIDE1_Non_EN_Discrete_AUC_test.png", height = 8, width = 8)

#============================================================================#
# plot predicted classes
#============================================================================#

ObsPred_CV_all <- data.frame(Score = ObsPred_CV_highRisk$ObservedScore,
                             Class = rep("Intermediate Risk (4-7)",nrow(ObsPred_CV_highRisk)))
ObsPred_CV_all$Class[(ObsPred_CV_highRisk$PredictedClass == "High") &
                       (ObsPred_CV_lowRisk$PredictedClass == "High/Intermediate")] <- "High Risk (8-14)"
ObsPred_CV_all$Class[(ObsPred_CV_lowRisk$PredictedClass == "Low") &
                       (ObsPred_CV_highRisk$PredictedClass == "Low/Intermediate")] <- "Low Risk (0-3)"

ObsPred_CV_all$Class <- factor(ObsPred_CV_all$Class,
                               levels = c("Low Risk (0-3)",
                                          "Intermediate Risk (4-7)",
                                          "High Risk (8-14)"))
p <- ggplot(ObsPred_CV_all) +
  geom_boxplot(aes(x = Class, y = Score, fill = Class), alpha = 0.5) +
  geom_jitter(aes(x= Class, y = Score, color = Class), width = 0.2, height = 0.1, alpha = 0.2) +
  xlab("Predicted Class") +
  ylab("CAIDE1 Score") +
  ggtitle("Performance in Cross-Validation") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, "Reds")[2:4]) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, "Reds")[2:4]) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 


ggsave(p, file = "CAIDE1_Non_EN_Discrete_Boxplot.png", height = 6, width = 8)


ObsPred_test_all <- data.frame(Score = ObsPred_test_highRisk$ObservedScore,
                             Class = rep("Intermediate Risk (4-7)",nrow(ObsPred_test_highRisk)))
ObsPred_test_all$Class[(ObsPred_test_highRisk$PredictedClass == "High") &
                       (ObsPred_test_lowRisk$PredictedClass == "High/Intermediate")] <- "High Risk (8-14)"
ObsPred_test_all$Class[(ObsPred_test_lowRisk$PredictedClass == "Low") &
                       (ObsPred_test_highRisk$PredictedClass == "Low/Intermediate")] <- "Low Risk (0-3)"

ObsPred_test_all$Class <- factor(ObsPred_test_all$Class,
                               levels = c("Low Risk (0-3)",
                                          "Intermediate Risk (4-7)",
                                          "High Risk (8-14)"))
p <- ggplot(ObsPred_test_all) +
  geom_boxplot(aes(x = Class, y = Score, fill = Class), alpha = 0.5) +
  geom_jitter(aes(x= Class, y = Score, color = Class), width = 0.2, height = 0.1, alpha = 1) +
  xlab("Predicted Class") +
  ylab("CAIDE1 Score") +
  ggtitle("Performance in Test Set") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, "Reds")[2:4]) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, "Reds")[2:4]) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 


ggsave(p, file = "CAIDE1_Non_EN_Discrete_Boxplot_test.png", height = 6, width = 8)



