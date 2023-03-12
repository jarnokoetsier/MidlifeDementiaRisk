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

# Selected features without prior feature filtering
Score = "CAIDE1"
methods = "Non"
load(paste0("CV_CAIDE1/FeatureSelection/CV_", Score, "_", methods,".RData"))

selectedFeatures <- rownames(as.matrix(coef(finalModel)))[as.matrix(coef(finalModel))[,1] != 0]
selectedFeatures <- selectedFeatures[-1]

# Correlation-based feature selection
FeatureSelection = "Cor_CAIDE1"
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/",f))
}

CorFeatures <- rownames(X_test_Cor)

length(intersect(CorFeatures, selectedFeatures))

# Selected features with prior correlation-based feature filtering
Score = "CAIDE1"
methods = "Cor"
load(paste0("CV_CAIDE1/CV_", Score, "_", methods,"_EN.RData"))

selectedFeatures_Cor <- rownames(varImp(finalModel, scale = FALSE)$importance)[varImp(finalModel, scale = FALSE)$importance != 0]

length(intersect(selectedFeatures_Cor, selectedFeatures))
length(intersect(CorFeatures, selectedFeatures_Cor))
