# Clear workspace and console
rm(list = ls())
cat("\014") 

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
###############################################################################

# All

###############################################################################

factors <- c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", "Physical", "HeartDisease",
             "Education", "Depression", "SysBP", "Diet")
factorNames <- c("BMI", "Type II Diabetes", "L-M Alchol Intake", "HDL Chol.",
                 "Total Chol.", "Physical Act.", "Heart Disease", "Education",
                 "Depression", "Systolic BP", "Dietary Intake")

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

colors <- c(RColorBrewer::brewer.pal(n = 6,"Dark2"),
            RColorBrewer::brewer.pal(n = 6,"Set2"))
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

ggsave(p, file = "AUC_all_Methylation.png", width = 8, height = 6)
  
###############################################################################

# BMI

###############################################################################

# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$BMI$AUCplot
auc_test_lit <- finalOutput_test$BMI$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$BMI$AUCplot
auc_CV_lit <- finalOutput$BMI$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$BMI$AUCplot
auc_test_EN <- finalOutput_test$BMI$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$BMI$AUCplot
auc_CV_EN <- finalOutput$BMI$AUC


# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$BMI$AUCplot
auc_test_RF <- finalOutput_test$BMI$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$BMI$AUCplot
auc_CV_RF <- finalOutput$BMI$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
               aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", round(auc_test_EN,2)), 
            color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", round(auc_test_RF,2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", round(auc_test_lit,2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("BMI", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Oranges")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/BMI_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", round(auc_CV_EN,2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", round(auc_CV_RF,2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", round(auc_CV_lit,2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("BMI", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Oranges")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/BMI_AUCplot_CV.png", width = 8, height = 6)



###############################################################################

# Diabetes

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$Diabetes$AUCplot
auc_test_lit <- finalOutput_test$Diabetes$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$Diabetes$AUCplot
auc_CV_lit <- finalOutput$Diabetes$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$Diabetes$AUCplot
auc_test_EN <- finalOutput_test$Diabetes$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$Diabetes$AUCplot
auc_CV_EN <- finalOutput$Diabetes$AUC


# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$Diabetes$AUCplot
auc_test_RF <- finalOutput_test$Diabetes$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$Diabetes$AUCplot
auc_CV_RF <- finalOutput$Diabetes$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Type II Diabetes", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "PuRd")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/T2D_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Type II Diabetes", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "PuRd")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/T2D_AUCplot_CV.png", width = 8, height = 6)



###############################################################################

# L-M Alcohol

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$Alcohol$AUCplot
auc_test_lit <- finalOutput_test$Alcohol$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$Alcohol$AUCplot
auc_CV_lit <- finalOutput$Alcohol$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$Alcohol$AUCplot
auc_test_EN <- finalOutput_test$Alcohol$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$Alcohol$AUCplot
auc_CV_EN <- finalOutput$Alcohol$AUC

# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$Alcohol$AUCplot
auc_test_RF <- finalOutput_test$Alcohol$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$Alcohol$AUCplot
auc_CV_RF <- finalOutput$Alcohol$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("L-M Alcohol Intake", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Blues")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Alcohol_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("L-M Alcohol Intake", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Blues")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Alcohol_AUCplot_CV.png", width = 8, height = 6)


###############################################################################

# HDL

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$HDL$AUCplot
auc_test_lit <- finalOutput_test$HDL$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$HDL$AUCplot
auc_CV_lit <- finalOutput$HDL$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$HDL$AUCplot
auc_test_EN <- finalOutput_test$HDL$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$HDL$AUCplot
auc_CV_EN <- finalOutput$HDL$AUC

# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$HDL$AUCplot
auc_test_RF <- finalOutput_test$HDL$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$HDL$AUCplot
auc_CV_RF <- finalOutput$HDL$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("HDL Cholesterol", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Reds")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/HDL_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("HDL Cholesterol", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Reds")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/HDL_AUCplot_CV.png", width = 8, height = 6)


###############################################################################

# Total Cholesterol

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$TotalChol$AUCplot
auc_test_lit <- finalOutput_test$TotalChol$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$TotalChol$AUCplot
auc_CV_lit <- finalOutput$TotalChol$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$TotalChol$AUCplot
auc_test_EN <- finalOutput_test$TotalChol$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$TotalChol$AUCplot
auc_CV_EN <- finalOutput$TotalChol$AUC

# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$TotalChol$AUCplot
auc_test_RF <- finalOutput_test$TotalChol$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$TotalChol$AUCplot
auc_CV_RF <- finalOutput$TotalChol$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Total Cholesterol", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "YlOrBr")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/TotalChol_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Total Cholesterol", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "YlOrBr")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/TotalChol_AUCplot_CV.png", width = 8, height = 6)


###############################################################################

# Physical Activity

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$Physical$AUCplot
auc_test_lit <- finalOutput_test$Physical$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$Physical$AUCplot
auc_CV_lit <- finalOutput$Physical$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$Physical$AUCplot
auc_test_EN <- finalOutput_test$Physical$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$Physical$AUCplot
auc_CV_EN <- finalOutput$Physical$AUC

# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$Physical$AUCplot
auc_test_RF <- finalOutput_test$Physical$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$Physical$AUCplot
auc_CV_RF <- finalOutput$Physical$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Physical Activity", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Greens")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Physical_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Physical Activity", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Greens")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Physical_AUCplot_CV.png", width = 8, height = 6)


###############################################################################

# Heart Disease

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$HeartDisease$AUCplot
auc_test_lit <- finalOutput_test$HeartDisease$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$HeartDisease$AUCplot
auc_CV_lit <- finalOutput$HeartDisease$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$HeartDisease$AUCplot
auc_test_EN <- finalOutput_test$HeartDisease$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$HeartDisease$AUCplot
auc_CV_EN <- finalOutput$HeartDisease$AUC

# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$HeartDisease$AUCplot
auc_test_RF <- finalOutput_test$HeartDisease$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$HeartDisease$AUCplot
auc_CV_RF <- finalOutput$HeartDisease$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Heart Disease", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "PuBuGn")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/HeartDisease_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Heart Disease", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "PuBuGn")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/HeartDisease_AUCplot_CV.png", width = 8, height = 6)


###############################################################################

# Education

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$Education$AUCplot
auc_test_lit <- finalOutput_test$Education$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$Education$AUCplot
auc_CV_lit <- finalOutput$Education$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$Education$AUCplot
auc_test_EN <- finalOutput_test$Education$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$Education$AUCplot
auc_CV_EN <- finalOutput$Education$AUC

# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$Education$AUCplot
auc_test_RF <- finalOutput_test$Education$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$Education$AUCplot
auc_CV_RF <- finalOutput$Education$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Education", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Purples")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Education_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Education", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Purples")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Education_AUCplot_CV.png", width = 8, height = 6)


###############################################################################

# Depression

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$Depression$AUCplot
auc_test_lit <- finalOutput_test$Depression$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$Depression$AUCplot
auc_CV_lit <- finalOutput$Depression$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$Depression$AUCplot
auc_test_EN <- finalOutput_test$Depression$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$Depression$AUCplot
auc_CV_EN <- finalOutput$Depression$AUC

# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$Depression$AUCplot
auc_test_RF <- finalOutput_test$Depression$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$Depression$AUCplot
auc_CV_RF <- finalOutput$Depression$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Depression", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "RdPu")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Depression_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Depression", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "RdPu")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Depression_AUCplot_CV.png", width = 8, height = 6)


###############################################################################

# Systolic blood pressure

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$SysBP$AUCplot
auc_test_lit <- finalOutput_test$SysBP$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$Depression$AUCplot
auc_CV_lit <- finalOutput$Depression$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$SysBP$AUCplot
auc_test_EN <- finalOutput_test$SysBP$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$SysBP$AUCplot
auc_CV_EN <- finalOutput$SysBP$AUC

# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$SysBP$AUCplot
auc_test_RF <- finalOutput_test$SysBP$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$SysBP$AUCplot
auc_CV_RF <- finalOutput$SysBP$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Systolic Blood Pressure", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Greys")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/SysBP_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Systolic Blood Pressure", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 6, "Greys")[4:6]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/SysBP_AUCplot_CV.png", width = 8, height = 6)


###############################################################################

# Diet

###############################################################################


# Literature-based feature selection + EN
load("~/PerFactor/Literature/finalOutput_lit_test.RData")
plotAUC_test_lit <- finalOutput_test$Diet$AUCplot
auc_test_lit <- finalOutput_test$Diet$AUC

load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
plotAUC_CV_lit <- finalOutput$Diet$AUCplot
auc_CV_lit <- finalOutput$Diet$AUC

# Correlation-based feature selection + EN
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
plotAUC_test_EN <- finalOutput_test$Diet$AUCplot
auc_test_EN <- finalOutput_test$Diet$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
plotAUC_CV_EN <- finalOutput$Diet$AUCplot
auc_CV_EN <- finalOutput$Diet$AUC

# Correlation-based feature selection + RF
load("~/PerFactor/Random Forest/finalOutput_test.RData")
plotAUC_test_RF <- finalOutput_test$Diet$AUCplot
auc_test_RF <- finalOutput_test$Diet$AUC

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
plotAUC_CV_RF <- finalOutput$Diet$AUCplot
auc_CV_RF <- finalOutput$Diet$AUC


p <- ggplot() +
  geom_path(data = plotAUC_test_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_test_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_test_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_test_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_test_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Dietary Intake", subtitle = "Performance in Test Set") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 11, "Spectral")[1:3]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Diet_AUCplot_test.png", width = 8, height = 6)


p <- ggplot() +
  geom_path(data = plotAUC_CV_EN, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, ElasticNet"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_RF, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Correlation, Random Forest"),
            linewidth = 1.5) +
  geom_path(data = plotAUC_CV_lit, 
            aes(y = Sensitivity, x = 1- Specificity, color = "Literature, ElasticNet"),
            linewidth = 1.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.85, y = 0.1, label = paste0("AUC: ", format(round(auc_CV_EN,2), nsmall = 2)), 
                color = "Correlation, ElasticNet"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.15, label = paste0("AUC: ", format(round(auc_CV_RF,2), nsmall = 2)), 
                color = "Correlation, Random Forest"),fontface = "bold") +
  geom_text(aes(x = 0.85, y = 0.2, label = paste0("AUC: ", format(round(auc_CV_lit,2), nsmall = 2)), 
                color = "Literature, ElasticNet"),fontface = "bold") +
  ggtitle("Dietary Intake", subtitle = "Performance in Cross-validation") +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 11, "Spectral")[1:3]) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, filename = "PerFactor/Diet_AUCplot_CV.png", width = 8, height = 6)