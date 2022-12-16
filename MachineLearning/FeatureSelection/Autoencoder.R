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
FeatureSelection = "S"
files <- list.files(paste0("X_", FeatureSelection))
for (f in files){
  load(paste0("X_", FeatureSelection, "/",f))
}


h2o.init()
features <- as.h2o(t(X_test_S))

ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = 100,
  activation = 'Tanh',
  sparse = FALSE
)

ae1_codings <- h2o.deepfeatures(ae1, features, layer = 1)
test <- as.data.frame(ae1_codings)


pred <- h2o.predict(ae1_codings, as.h2o(X_nonTest[,1:100]))
