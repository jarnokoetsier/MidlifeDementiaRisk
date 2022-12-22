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
FeatureSelection = "var"

# Load CV index
load("CVindex_CAIDE1.RData")


################################################################################

# ROC in CV

################################################################################


#****************************************************************************#
# Cross-validation
#****************************************************************************#

#============================================================================#
# Low risk
#============================================================================#
# Load low risk model training
load(paste0("CV_", Score, "_", FeatureSelection,"_LowRisk.RData"))

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
#============================================================================#
# High risk
#============================================================================#

# Load high risk model
load(paste0("CV_", Score, "_", FeatureSelection,"_HighRisk.RData"))

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

#============================================================================#
# plot ROC
#============================================================================#

ggplot() +
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
ggplot(ObsPred_CV_all) +
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




#****************************************************************************#
# Test set
#****************************************************************************#


# Get performance on test set
X_test_var_M <- log2(X_test_var/(1-X_test_var))
test_prob <- predict(finalModel, t(X_test_var_M), type = "response")[,1]
test_class <-factor(ifelse(Y_test$CAIDE <4, "Low", "Intermediate_High"),
                               levels = c("Intermediate_High", "Low"))

roc_list_test <- roc(response = test_class, 
                     predictor = test_prob)

plot(roc_list_test)
auc(roc_list_test)


output[["High"]] <- roc_list_test



# Load data
files <- list.files(paste0("X_", FeatureSelection))
for (f in files){
  load(paste0("X_", FeatureSelection, "/", f))
}

# Prepare data
X_train = log2(X_CAIDE1_var/(1-X_CAIDE1_var))
Y_train = Y_CAIDE1$CAIDE

# Test if samples are in correct order
all(colnames(X_train) == Y_CAIDE1$Basename)


# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "var"
