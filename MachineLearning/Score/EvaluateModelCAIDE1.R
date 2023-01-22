# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load cell type composition and phenotype data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")


################################################################################

# Compare performances of feature selection approaches

################################################################################

# Score and feature selection method
Score = "CAIDE1"
methods = c("S", "var", "varM", "varCor", "varMCor","PC", "KS", "Non")

# Load data
for (i in 1:length(methods)){
  load(paste0("X/X_", methods[i], "/", "X_test_", methods[i], ".RData"))
}
# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Retrieve performance in CV for the different feature selection methods
Performances <- list()
Performances_test <- list()
for (i in 1:length(methods)){
  
  # load output
  load(paste0("CV_CAIDE1/CV_", Score, "_", methods[i],".RData"))
  
  # Put performance into list
  Performances[[i]] <- perf
  
  # Performance in test data
  X_test <- get(paste0("X_test_", methods[i]))
  if (methods[i] != "PC"){
    testData <- log2(X_test/(1-X_test))
  }
  if (methods[i] == "PC"){
    testData <- t(X_test)
  }
  pred_test <- predict(finalModel, t(testData))
  Performances_test[[i]] <- RMSE(pred = pred_test,
                                 obs = Y_test$CAIDE)
}
names(Performances) <- methods
names(Performances_test) <- methods



# Combine into single data frame
performanceDF <- data.frame(
  RMSE = unlist(Performances),
  FeatureSelection = c(rep("S-score", 25),
                       rep("Variance (\u03b2)", 25),
                       rep("Variance (M)", 25),
                       rep("Variance (\u03b2, Cor)",25),
                       rep("Variance (M, Cor)", 25),
                       rep("PCA",25),
                       rep("KS-like",25),
                       rep("None",25))
)
performanceDF$FeatureSelection <- factor(performanceDF$FeatureSelection,
                                         levels = unique(performanceDF$FeatureSelection))

# Combine into single data frame
performanceDF_test <- data.frame(
  RMSE = unlist(Performances_test),
  FeatureSelection = c(rep("S-score", 1),
                       rep("Variance (\u03b2)", 1),
                       rep("Variance (M)", 1),
                       rep("Variance (\u03b2, Cor)",1),
                       rep("Variance (M, Cor)", 1),
                       rep("PCA",1),
                       rep("KS-like",1),
                       rep("None", 1))
)
performanceDF_test$FeatureSelection <- factor(performanceDF_test$FeatureSelection,
                                              levels = unique(performanceDF_test$FeatureSelection))



# Load data
load(paste0("CV_CAIDE1/CV_CAIDE1_Cor_EN.RData"))
load("X/X_Cor_CAIDE1/X_test_Cor.RData")

# Performance in training data
optPar <- which.min(rowMeans(perf))
optPerf <- NULL
for (i in 1:length(trainResults)){
  optPerf <- c(optPerf,trainResults[[i]]$RMSE[optPar])
}
performanceDF <- rbind.data.frame(performanceDF, data.frame(RMSE = optPerf,
                                                            FeatureSelection = rep("Correlation", 25)))

# Performance in test data
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred_test <- predict(finalModel, t(testData))
perf_test <- RMSE(pred = pred_test,
                  obs = Y_test$CAIDE)

performanceDF_test <- rbind.data.frame(performanceDF_test, data.frame(RMSE = perf_test,
                                                                      FeatureSelection = rep("Correlation", 1)))


performanceDF$FeatureSelection <- factor(performanceDF$FeatureSelection,
                                         levels = c(rep("S-score", 1),
                                                    rep("Variance (\u03b2)", 1),
                                                    rep("Variance (M)", 1),
                                                    rep("Variance (\u03b2, Cor)",1),
                                                    rep("Variance (M, Cor)", 1),
                                                    rep("PCA",1),
                                                    rep("KS-like",1),
                                                    rep("Correlation",1),
                                                    rep("None", 1)))


# Make plot
p <- ggplot(performanceDF) +
  geom_boxplot(aes(x = FeatureSelection, y = RMSE, fill = FeatureSelection), alpha = 0.3) +
  geom_point(aes(x = FeatureSelection, y = RMSE, color = FeatureSelection), 
             position=position_jitterdodge(jitter.width = 1), size = 2) +
  geom_point(data = performanceDF_test, aes(x = FeatureSelection, y = RMSE), 
             color = "black", size = 5, shape = 18, alpha = 0.7) +
  xlab("Feature Selection Method") +
  ggtitle(Score) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2")[c(1:5,8,6,7)], "red")) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2")[c(1:5,8,6,7)], "red")) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = paste0(Score, "_RMSE_boxplot.png"), width = 10, height = 6)


