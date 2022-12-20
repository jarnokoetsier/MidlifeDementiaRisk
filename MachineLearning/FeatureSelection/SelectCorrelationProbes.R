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
source("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA/MachineLearning/FUN_MachineLearning.R")

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load cell type composition and phenotype data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")


# Meta data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

load("E:/Thesis/EXTEND/Methylation/FeatureSelection/correlations_CAIDE1factors_named.RData")

load(paste0("X_cor/","X_cor.RData"))

X_nonTest_cor <- X_cor[,Y_CAIDE1$Basename]

all(colnames(X_nonTest_cor) == Y_CAIDE1$Basename)
factors <- Y_CAIDE1[,14:21]

correlations1 <- matrix(NA, nrow = nrow(X_nonTest_cor), ncol = ncol(factors))
for (i in 1:ncol(factors)) {
  correlations1[,i] <- apply(X_nonTest_cor, 1, 
                            function(x){cor(x, 
                                            factors[,i], 
                                            method = "spearman")})
}

plot(correlations1[,1], correlations[rownames(X_nonTest_cor),1])



X_nonTest_cor <- X_cor[,Y_nonTest$Basename]
X_CAIDE1_cor <- X_cor[,Y_CAIDE1$Basename]
X_LIBRA_cor <- X_cor[,Y_LIBRA$Basename]
X_test_cor <- X_cor[,Y_test$Basename]

save(X_nonTest_cor, file = "X_nonTest_cor.RData")
save(X_CAIDE1_cor, file = "X_CAIDE1_cor.RData")
save(X_LIBRA_cor, file = "X_LIBRA_cor.RData")
save(X_test_cor, file = "X_test_cor.RData")



# Select top correlated features for each factor
probes <- list()
for (i in 1:7){
  probes[[i]] <- names(tail(sort(abs(correlations[,i])),1700))
}

# get exactly 10,000 probes
n = 1
output <- unique(unlist(probes))
while (length(output) > 10000){
  probes[[n]] <- probes[[n]][-1]
  output <- unique(unlist(probes))
  
  if (n < 7){
    n = n + 1
  } else {
    n = 1
  }
}





# Feature selection data
FeatureSelection = "S"
files <- list.files(paste0("X_", FeatureSelection))
for (f in files){
  load(paste0("X_", FeatureSelection, "/",f))
}

