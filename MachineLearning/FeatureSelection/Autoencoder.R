library(h2o)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Cell type composition data
load("cellType.RData")

# Phenotype data
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Meta data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Feature selection data
FeatureSelection = "PC"
files <- list.files(paste0("X_", FeatureSelection))
for (f in files){
  load(paste0("X_", FeatureSelection, "/",f))
}


h2o.init()
features <- as.h2o(X_nonTest_PC)

ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = c(200,50,200),
  activation = 'Tanh',
  sparse = FALSE
)

ae1_codings <- h2o.deepfeatures(ae1, features, layer = 1)
test <- as.data.frame(ae1_codings)


pred <- h2o.predict(ae1, as.h2o(X_CAIDE1_PC))
pred <- as.data.frame(pred)


hyper_grid <- list(hidden = list(
  c(50),
  c(100), 
  c(300, 100, 300),
  c(100, 50, 100),
  c(250, 100, 50, 100, 250)
))

# Execute grid search
ae_grid <- h2o.grid(
  algorithm = 'deeplearning',
  x = seq_along(features),
  training_frame = features,
  grid_id = 'autoencoder_grid',
  autoencoder = TRUE,
  activation = 'Tanh',
  hyper_params = hyper_grid,
  sparse = TRUE,
  ignore_const_cols = FALSE,
  seed = 123
)


grid <- h2o.getGrid('autoencoder_grid', sort_by = 'mse', decreasing = FALSE)
best_model_id <- grid@model_ids[[5]]
best_model <- h2o.getModel(best_model_id)
pred <- h2o.predict(best_model, as.h2o(X_CAIDE1_PC))
pred <- as.data.frame(pred)
