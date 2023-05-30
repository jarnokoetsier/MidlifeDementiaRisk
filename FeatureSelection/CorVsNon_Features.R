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

# Selected features without prior feature filtering:

# 1. Load EN model w/o feature selection
Score = "CAIDE1"
methods = "Non"
load(paste0("CV_CAIDE1/FeatureSelection/CV_", Score, "_", methods,".RData"))

# 2. Retrieve features from model
selectedFeatures <- rownames(as.matrix(coef(finalModel)))[as.matrix(coef(finalModel))[,1] != 0]
selectedFeatures <- selectedFeatures[-1]

# Selected features by correlation-based feature selection

# 1. load features selected by correlation-based feature selection
FeatureSelection = "Cor_CAIDE1"
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/",f))
}

# 2. collect selected features
CorFeatures <- rownames(X_test_Cor)

# Selected features with prior correlation-based feature filtering

# 1. Load EN model with feature selection
Score = "CAIDE1"
methods = "Cor"

# 2. Retrieve features from model
load(paste0("CV_CAIDE1/CV_", Score, "_", methods,"_EN.RData"))
selectedFeatures_Cor <- rownames(varImp(finalModel, scale = FALSE)$importance)[varImp(finalModel, scale = FALSE)$importance != 0]

# Compare common features
length(intersect(selectedFeatures_Cor, selectedFeatures))
length(intersect(CorFeatures, selectedFeatures_Cor))
length(intersect(CorFeatures, selectedFeatures))
