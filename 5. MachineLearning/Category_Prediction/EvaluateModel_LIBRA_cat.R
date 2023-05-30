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
files <- list.files('Y', pattern = ".RData")
for (f in files){
  load(paste0("Y/",f))
}


################################################################################

# ROC

################################################################################

methods <- c("EN", "sPLSDA", "RF")
methodNames <- c("ElasticNet", "sPLS-DA", "Random Forest")

# Low Risk
plotDF_LowRisk_CV <- NULL
plotDF_LowRisk_test <- NULL
aucDF_LowRisk_CV <- NULL
aucDF_LowRisk_test <- NULL
for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_LIBRA/Performance_LIBRA_LowRisk_Cor_", methods[i],".RData"))
  
  # construct ROC
  roc_list_CV <- roc(response = factor(ObsPred_CV_LowRisk$obs,
                                       levels = c("Intermediate_High", "Low")), 
                     predictor = ObsPred_CV_LowRisk$p)
  
  # Get sensitivies and specificities
  temp <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                     Specificity = roc_list_CV$specificities,
                     Method = paste0(rep(methodNames[i], length(roc_list_CV$specificities))," (Low Risk)"))
  
  plotDF_LowRisk_CV <- rbind.data.frame(plotDF_LowRisk_CV, temp)
  
  temp_auc <- data.frame(AUC = as.numeric(auc(roc_list_CV)),
                         Method = paste0(methodNames[i]," (Low Risk)")
  )
  
  aucDF_LowRisk_CV <- rbind.data.frame(aucDF_LowRisk_CV, temp_auc)
  
  
  roc_list_test <- roc(response = factor(ObsPred_test_LowRisk$obs,
                                         levels = c("Intermediate_High", "Low")), 
                       predictor = ObsPred_test_LowRisk$p)
  
  temp <- data.frame(Sensitivity = roc_list_test$sensitivities,
                     Specificity = roc_list_test$specificities,
                     Method = paste0(rep(methodNames[i], length(roc_list_test$specificities))," (Low Risk)"))
  
  plotDF_LowRisk_test <- rbind.data.frame(plotDF_LowRisk_test, temp)
  
  temp_auc <- data.frame(AUC = as.numeric(auc(roc_list_test)),
                         Method = paste0(methodNames[i]," (Low Risk)")
  )
  
  aucDF_LowRisk_test <- rbind.data.frame(aucDF_LowRisk_test, temp_auc)
  
}

# High Risk
plotDF_HighRisk_CV <- NULL
plotDF_HighRisk_test <- NULL
aucDF_HighRisk_CV <- NULL
aucDF_HighRisk_test <- NULL
for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_LIBRA/Performance_LIBRA_HighRisk_Cor_", methods[i],".RData"))
  
  # construct ROC
  roc_list_CV <- roc(response = factor(ObsPred_CV_HighRisk$obs,
                                       levels = c("Low_Intermediate", "High")), 
                     predictor = ObsPred_CV_HighRisk$p)
  
  # Get sensitivies and specificities
  temp <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                     Specificity = roc_list_CV$specificities,
                     Method = paste0(rep(methodNames[i], length(roc_list_CV$specificities))," (High Risk)"))
  
  plotDF_HighRisk_CV <- rbind.data.frame(plotDF_HighRisk_CV, temp)
  
  temp_auc <- data.frame(AUC = as.numeric(auc(roc_list_CV)),
                         Method = paste0(methodNames[i]," (High Risk)")
  )
  
  aucDF_HighRisk_CV <- rbind.data.frame(aucDF_HighRisk_CV, temp_auc)
  
  
  roc_list_test <- roc(response = factor(ObsPred_test_HighRisk$obs,
                                         levels = c("Low_Intermediate", "High")), 
                       predictor = ObsPred_test_HighRisk$p)
  
  temp <- data.frame(Sensitivity = roc_list_test$sensitivities,
                     Specificity = roc_list_test$specificities,
                     Method = paste0(rep(methodNames[i], length(roc_list_test$specificities))," (High Risk)"))
  
  plotDF_HighRisk_test <- rbind.data.frame(plotDF_HighRisk_test, temp)
  
  temp_auc <- data.frame(AUC = as.numeric(auc(roc_list_test)),
                         Method = paste0(methodNames[i]," (High Risk)")
  )
  
  aucDF_HighRisk_test <- rbind.data.frame(aucDF_HighRisk_test, temp_auc)
  
}
plotDF_all_CV <- rbind.data.frame(plotDF_LowRisk_CV, plotDF_HighRisk_CV)
plotDF_all_CV$Method <- factor(plotDF_all_CV$Method,
                               levels = c("ElasticNet (Low Risk)",
                                          "sPLS-DA (Low Risk)",
                                          "Random Forest (Low Risk)",
                                          "ElasticNet (High Risk)",
                                          "sPLS-DA (High Risk)",
                                          "Random Forest (High Risk)"))
aucDF_all_CV <- rbind.data.frame(aucDF_LowRisk_CV, aucDF_HighRisk_CV)
aucDF_all_CV$Method <- factor(aucDF_all_CV$Method,
                              levels = c("ElasticNet (Low Risk)",
                                         "sPLS-DA (Low Risk)",
                                         "Random Forest (Low Risk)",
                                         "ElasticNet (High Risk)",
                                         "sPLS-DA (High Risk)",
                                         "Random Forest (High Risk)"))
plotDF_all_test <- rbind.data.frame(plotDF_LowRisk_test, plotDF_HighRisk_test)
plotDF_all_test$Method <- factor(plotDF_all_test$Method,
                                 levels = c("ElasticNet (Low Risk)",
                                            "sPLS-DA (Low Risk)",
                                            "Random Forest (Low Risk)",
                                            "ElasticNet (High Risk)",
                                            "sPLS-DA (High Risk)",
                                            "Random Forest (High Risk)"))
aucDF_all_test <- rbind.data.frame(aucDF_LowRisk_test, aucDF_HighRisk_test)
aucDF_all_test$Method <- factor(aucDF_all_test$Method,
                                levels = c("ElasticNet (Low Risk)",
                                           "sPLS-DA (Low Risk)",
                                           "Random Forest (Low Risk)",
                                           "ElasticNet (High Risk)",
                                           "sPLS-DA (High Risk)",
                                           "Random Forest (High Risk)"))

# Set colors
colors <- c(brewer.pal(n = 5, name = "Reds")[c(3)],
            brewer.pal(n = 5, name = "Reds")[c(4)],
            brewer.pal(n = 5, name = "Reds")[c(5)],
            brewer.pal(n = 5, name = "Blues")[c(3)],
            brewer.pal(n = 5, name = "Blues")[c(4)],
            brewer.pal(n = 5, name = "Blues")[c(5)]
            )

# Plot performance in cross-validation
p <- ggplot() +
  geom_path(data = plotDF_all_CV, aes(y = Sensitivity, x = 1- Specificity,
                                      color = Method), 
            linewidth = 1.5, linetype = "solid") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_text(data = aucDF_all_CV, aes(x = 0.9, y = rev(c(0.05,0.1,0.15,0.2,0.25, 0.3)),
                                     label = paste0("AUC: ", format(round(AUC,2), nsmall = 2)), 
                                     color = Method), fontface = "bold") +
  scale_color_manual(values = colors) +
  ggtitle("LIBRA") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "LIBRA_Cat_Cor_AUC_CV.png", height = 5, width = 7.5)

# Plot performance in test set
p <- ggplot() +
  geom_path(data = plotDF_all_test, aes(y = Sensitivity, x = 1- Specificity,
                                        color = Method), 
            linewidth = 1.5, linetype = "solid") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_text(data = aucDF_all_test, aes(x = 0.9, y = rev(c(0.05,0.1,0.15,0.2,0.25, 0.3)),
                                     label = paste0("AUC: ", format(round(AUC,2), nsmall = 2)), 
                                     color = Method), fontface = "bold") +
  scale_color_manual(values = colors) +
  ggtitle("LIBRA") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "LIBRA_Cat_Cor_AUC_test.png", height = 5, width = 7.5)


################################################################################

# Boxplots

################################################################################
methods <- c("EN", "sPLSDA", "RF")
methodNames <- c("ElasticNet", "sPLS-DA", "Random Forest")

ObsPred_CV_all <- NULL
ObsPred_test_all <- NULL
for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_LIBRA/Performance_LIBRA_HighRisk_Cor_", methods[i],".RData"))
  load(paste0("CV_LIBRA/Performance_LIBRA_LowRisk_Cor_", methods[i],".RData"))
  
  # Obs vs predicted in CV
  temp <- data.frame(Score = ObsPred_CV_HighRisk$ObservedScore,
                     Class = rep("Intermediate Risk (0-2)",nrow(ObsPred_CV_HighRisk)))
  temp$Class[(ObsPred_CV_HighRisk$pred == "High") &
               (ObsPred_CV_LowRisk$pred == "Intermediate_High")] <- "High Risk (> 2)"
  temp$Class[(ObsPred_CV_LowRisk$pred == "Low") &
               (ObsPred_CV_HighRisk$pred == "Low_Intermediate")] <- "Low Risk (< 0)"
  temp$Method <- rep(methodNames[i], nrow(temp))
  
  ObsPred_CV_all <- rbind.data.frame(ObsPred_CV_all, temp)
  
  # Obs vs predicted in test set
  temp <- data.frame(Score = ObsPred_test_HighRisk$ObservedScore,
                     Class = rep("Intermediate Risk (0-2)",nrow(ObsPred_test_HighRisk)))
  temp$Class[(ObsPred_test_HighRisk$pred == "High") &
               (ObsPred_test_LowRisk$pred == "Intermediate_High")] <- "High Risk (> 2)"
  temp$Class[(ObsPred_test_LowRisk$pred == "Low") &
               (ObsPred_test_HighRisk$pred == "Low_Intermediate")] <- "Low Risk (< 0)"
  temp$Method <- rep(methodNames[i], nrow(temp))
  
  ObsPred_test_all <- rbind.data.frame(ObsPred_test_all, temp)
  
}


ObsPred_CV_all$Class <- factor(ObsPred_CV_all$Class,
                               levels = c("Low Risk (< 0)","Intermediate Risk (0-2)","High Risk (> 2)"))
ObsPred_CV_all$Method <- factor(ObsPred_CV_all$Method,
                                levels = methodNames)

# Observed vs predicted in cross-validation
p <- ggplot(ObsPred_CV_all) +
  geom_rect(ymin = -Inf, ymax = 0, xmin = -Inf, xmax = Inf,fill = "#FFF5F0", alpha = 0.5) +
  geom_rect(ymin = 0, ymax = 2, xmin = -Inf, xmax = Inf,fill = "#FEE0D2", alpha = 0.5) +
  geom_rect(ymin = 2, ymax = Inf, xmin = -Inf, xmax = Inf,fill = "#FCBBA1", alpha = 0.5) +
  geom_hline(yintercept = 0, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 2, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_boxplot(aes(x = Class, y = Score, fill = Method), 
               alpha = 1,outlier.shape = NA) +
  geom_point(aes(x= Class, y = Score, color = Method, group = Method),
             position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1)) +
  xlab("Predicted Class") +
  ylab("LIBRA Score") +
  scale_y_continuous(breaks = c(-2,0,2,4,6,8,10,12), labels = c(-2,0,2,4,6,8,10,12))+
  ggtitle("LIBRA") +
  scale_fill_manual(values = c("#FB6A4A","#FEC44F", "#BCBDDC")) +
  scale_color_manual(values = c("#A50F15","#CC4C02", "#6A51A3")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

# Save plot
ggsave(p, file = "LIBRA_Cat_Cor_boxPlot_CV.png", height = 5, width = 7.5)



ObsPred_test_all$Class <- factor(ObsPred_test_all$Class,
                                 levels = c("Low Risk (< 0)","Intermediate Risk (0-2)","High Risk (> 2)"))
ObsPred_test_all$Method <- factor(ObsPred_test_all$Method,
                                  levels = methodNames)

# Observed vs predicted in test set
p <- ggplot(ObsPred_test_all) +
  geom_rect(ymin = -Inf, ymax = 0, xmin = -Inf, xmax = Inf,fill = "#FFF5F0", alpha = 0.5) +
  geom_rect(ymin = 0, ymax = 2, xmin = -Inf, xmax = Inf,fill = "#FEE0D2", alpha = 0.5) +
  geom_rect(ymin = 2, ymax = Inf, xmin = -Inf, xmax = Inf,fill = "#FCBBA1", alpha = 0.5) +
  geom_hline(yintercept = 0, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 2, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_boxplot(aes(x = Class, y = Score, fill = Method), 
               alpha = 1,outlier.shape = NA) +
  geom_point(aes(x= Class, y = Score, color = Method, group = Method),
             position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1)) +
  xlab("Predicted Class") +
  ylab("LIBRA Score") +
  ggtitle("LIBRA") +
  scale_y_continuous(breaks = c(-2,0,2,4,6,8,10), labels = c(-2,0,2,4,6,8,10))+
  scale_fill_manual(values = c("#FB6A4A","#FEC44F", "#BCBDDC")) +
  scale_color_manual(values = c("#A50F15","#CC4C02", "#6A51A3")) +
  #scale_fill_manual(values = c("#FB6A4A","#6BAED6", "#FEC44F", "#BCBDDC")) +
  #scale_color_manual(values = c("#A50F15","#08519C", "#CC4C02", "#6A51A3")) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

# Save plot
ggsave(p, file = "LIBRA_Cat_Cor_boxPlot_test.png", height = 5, width = 7.5)

# Make confusion matrix
model <- "sPLS-DA"
test <- ObsPred_test_all[ObsPred_test_all$Method == model,]
test$obs <- rep("Intermediate Risk (0-2)", nrow(test))
test$obs[test$Score < 0] <- "Low Risk (< 0)"
test$obs[test$Score > 2] <- "High Risk (> 2)"
confusionMatrix(factor(test$obs, levels = levels(test$Class)), test$Class)

################################################################################

# Which factors are being predicted

################################################################################

# Get factors of test set
load("X/X_Cor_LIBRA/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))

load("E:/Thesis/EXTEND/Phenotypes/EPILIBRA.Rdata")
rownames(EPILIBRA) <- EPILIBRA$Basename
Y_test_LIBRA <- EPILIBRA[colnames(testData),]


# tested methods
methods <- c("EN", "sPLSDA", "RF")
methodNames <- c("ElasticNet",  "sPLS-DA", "Random Forest")

# Low Risk
plotDF_all <- NULL
for (i in 1:length(methods)){
  # Load data
  load(paste0("CV_LIBRA/Performance_LIBRA_LowRisk_Cor_", methods[i],".RData"))
  load(paste0("CV_LIBRA/Performance_LIBRA_HighRisk_Cor_", methods[i],".RData"))
  
  # Combine observed and predicted
  dataMatrix <- cbind.data.frame(Y_test_LIBRA[,c(10:20)],factor(ObsPred_test_LowRisk$pred,
                                                                 levels = c("Intermediate_High", "Low")))
  dataMatrix <- cbind.data.frame(dataMatrix,factor(ObsPred_test_HighRisk$pred,
                                                   levels = c("Low_Intermediate", "High")))
  colnames(dataMatrix) <- c(colnames(Y_test_LIBRA[,c(10:20)]), "LowRisk", "HighRisk")
  
  Rsquared_low <- rep(NA, 11)
  Rsquared_high <- rep(NA, 11)
  Sig_low <- rep(NA, 11)
  Sig_high <- rep(NA, 11)
  Y_factors <- Y_test_LIBRA[,c(10:20)]
  for (factor in 1:ncol(Y_factors)){
    
    # Explained variance predicted score
    formula <- paste0("LowRisk ~ ", colnames(Y_factors)[factor])
    formula_null <- paste0("LowRisk ~ ",1)
    
    # Fit model
    model <- glm(as.formula(formula), data = as.data.frame(dataMatrix), family = "binomial")
    nullmodel <- glm(as.formula(formula_null), data = as.data.frame(dataMatrix), family = "binomial")
    
    # Calculate R-squared
    Rsquared_low[factor] = as.numeric(1-logLik(model)/(logLik(nullmodel)))
    Sig_low[factor] <- summary(model)$coefficients[2,4]
    
    #=========================================================================#
    
    # Explained variance observed score
    formula <- paste0("HighRisk ~ ", colnames(Y_factors)[factor])
    formula_null <- paste0("HighRisk ~ ",1)
    
    # Fit model
    model <- glm(as.formula(formula), data = as.data.frame(dataMatrix), family = "binomial")
    nullmodel <- glm(as.formula(formula_null), data = as.data.frame(dataMatrix), family = "binomial")
    
    # Calculate R-squared
    Rsquared_high[factor] = as.numeric(1-logLik(model)/(logLik(nullmodel)))
    Sig_high[factor] <- summary(model)$coefficients[2,4]
  }
  
  factorNames <- c("Healthy Diet", "Physical Inactivity", "Smoking", "L-M Alcohol Intake",
                   "BMI", "Depression", "Type 2 Diabetes","Systolic Blood Pressure", 
                   "HDL Cholesterol", "Heart Disease", "Kidney Disease")
  plotDF <- data.frame(Name = rep(factorNames,2), 
                       R2  = c(Rsquared_low, Rsquared_high),
                       Sig = c(Sig_low, Sig_high),
                       Model = c(rep("Low Risk Model",length(factorNames)),
                                 rep("High Risk Model",length(factorNames))),
                       Method = rep(methodNames[i],2*length(factorNames)))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

plotDF_all$Model <- factor(plotDF_all$Model,levels = c("Low Risk Model", "High Risk Model"))
plotDF_all$Method <- factor(plotDF_all$Method,
                            levels = methodNames)
plotDF_all$X <- plotDF_all$R2 + ifelse(plotDF_all$Sig < 0.05, 0.01,NA)

# Make plot
p <- ggplot(plotDF_all)+
  geom_bar(aes(x = Name, y = R2, fill = Method),
           stat = "identity", position=position_dodge(), color = "black") +
  geom_point(aes(x = Name, y = X, color = Method),
             position=position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.9),
             shape = 18) +
  coord_flip() +
  facet_grid(cols = vars(Model)) +
  xlab(NULL) +
  ylab(expression("McFadden's "~R^2)) +
  ylim(c(0,0.5)) +
  ggtitle("LIBRA") +
  labs(fill = NULL) +
  guides(color = "none") +
  theme_minimal() +
  scale_fill_manual(values = c("#EF3B2C","#FE9929", "#807DBA")) +
  scale_color_manual(values = c("black","black", "black")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

# Save plot
ggsave(p, file = "LIBRA_Cat_Cor_whichFactors_test.png", height = 6, width = 10)
