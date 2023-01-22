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

# Set direcotries
DataDir <- "Data/"
OutputDir <- "PreProcessing/"

# Load machine learning functions
source("FUN_MachineLearning.R")


################################################################################

# Correlation-based feature selection

################################################################################


#*****************************************************************************#
# CAIDE1
#*****************************************************************************#

# Load data
load(paste0(DataDir,"X_CAIDE1.Rdata"))
load(paste0(DataDir,"Y_CAIDE1.Rdata"))
load(paste0(DataDir,"X_nonTest.RData"))
load(paste0(DataDir,"Y_nonTest.RData"))
load(paste0(DataDir,"X_test.RData"))
load(paste0(DataDir,"Y_test.RData"))

# Check whether samples are in correct order
all(colnames(X_CAIDE1) == Y_CAIDE1$Basename)

# CAIDE1 factors
factors <- Y_CAIDE1[,14:20]

# Calculate spearman correlation
correlations <- matrix(NA, nrow = nrow(X_CAIDE1), ncol = ncol(factors))
for (i in 1:ncol(factors)) {
  correlations[,i] <- apply(X_CAIDE1, 1, 
                            function(x){cor(x, 
                                            factors[,i], 
                                            method = "spearman")})
}
rownames(correlations) <- rownames(X_CAIDE1)
colnames(correlations) <- colnames(factors)
save(correlations, file = paste0(DataDir,"correlations_CAIDE1factors.RData"))

# Select top correlated features for each factor
probes <- list()
for (p in 1:7){
  probes[[p]] <- names(tail(sort(abs(correlations[,p])),1700))
}

# get exactly 10,000 probes
n = 1
finalProbes <- unique(unlist(probes))
while (length(finalProbes) > 10000){
  probes[[n]] <- probes[[n]][-1]
  finalProbes <- unique(unlist(probes))
  
  if (n < 7){
    n = n + 1
  } else {
    n = 1
  }
}

X_nonTest_Cor <- X_nonTest[finalProbes, ]
X_CAIDE1_Cor <- X_CAIDE1[finalProbes, ]
X_test_Cor <- X_test[finalProbes, ]

save(X_nonTest_Cor, file = "CAIDE1_Cor/X_nonTest_Cor.RData")
save(X_CAIDE1_Cor, file = "CAIDE1_Cor/X_CAIDE1_Cor.RData")
save(X_test_Cor, file = "CAIDE1_Cor/X_test_Cor.RData")

# Probes for CV
CVProbes <- list()
for (j in colnames(correlations)){
  CVProbes[[j]] <- names(tail(sort(abs(correlations[,j])),10000))
}
CVProbes <- unique(unlist(CVProbes))

X_CAIDE1_CorCV <- X_CAIDE1[CVProbes, ]
save(X_CAIDE1_CorCV, file = "CAIDE1_Cor/X_CAIDE1_CorCV.RData")


#*****************************************************************************#
# LIBRA
#*****************************************************************************#

# Load data
load(paste0(DataDir,"X_LIBRA.RData"))
load(paste0(DataDir,"Y_LIBRA.RData"))
load(paste0(DataDir,"X_nonTest.RData"))
load(paste0(DataDir,"Y_nonTest.RData"))
load(paste0(DataDir,"X_test.RData"))
load(paste0(DataDir,"Y_test.RData"))

# Check whether samples are in correct order
all(colnames(X_LIBRA) == Y_LIBRA$Basename)

# CAIDE1 factors
factors <- Y_LIBRA[,14:24]

# Calculate spearman correlation
correlations <- matrix(NA, nrow = nrow(X_LIBRA), ncol = ncol(factors))
for (i in 1:ncol(factors)) {
  correlations[,i] <- apply(X_LIBRA, 1, 
                            function(x){cor(x, 
                                            factors[,i], 
                                            method = "spearman")})
}
rownames(correlations) <- rownames(X_LIBRA)
colnames(correlations) <- colnames(factors)
save(correlations, file = "LIBRA_Cor/correlations_LIBRAfactors.RData")

# Select top correlated features for each factor
probes <- list()
for (p in 1:7){
  probes[[p]] <- names(tail(sort(abs(correlations[,p])),1050))
}

# get exactly 10,000 probes
n = 1
finalProbes <- unique(unlist(probes))
while (length(finalProbes) > 10000){
  probes[[n]] <- probes[[n]][-1]
  finalProbes <- unique(unlist(probes))
  
  if (n < 7){
    n = n + 1
  } else {
    n = 1
  }
}

# Extract data
X_nonTest_CorL <- X_nonTest[finalProbes, ]
X_LIBRA_CorL <- X_CAIDE1[finalProbes, ]
X_test_CorL <- X_test[finalProbes, ]

# Save data
save(X_nonTest_CorL, file = "LIBRA_Cor/X_nonTest_CorL.RData")
save(X_LIBRA_CorL, file = "LIBRA_Cor/X_LIBRA_CorL.RData")
save(X_test_CorL, file = "LIBRA_Cor/X_test_CorL.RData")

# Probes for CV
CVProbes <- list()
for (j in colnames(correlations)){
  CVProbes[[j]] <- names(tail(sort(abs(correlations[,j])),8000))
}
CVProbes <- unique(unlist(CVProbes))

X_LIBRA_CorLCV <- X_LIBRA[CVProbes, ]
save(X_LIBRA_CorCV, file = "LIBRA_Cor/X_LIBRA_CorLCV.RData")


#*****************************************************************************#
# CAIDE1
#*****************************************************************************#

# Load data
load(paste0(DataDir,"X_CAIDE2.Rdata"))
load(paste0(DataDir,"Y_CAIDE2.RData"))
load(paste0(DataDir,"X_nonTest.RData"))
load(paste0(DataDir,"Y_nonTest.RData"))
load(paste0(DataDir,"X_test.RData"))
load(paste0(DataDir,"Y_test.RData"))

# Check whether samples are in correct order
all(colnames(X_CAIDE2) == Y_CAIDE2$Basename)

# CAIDE1 factors
factors <- Y_CAIDE2[,14:21]

# Calculate spearman correlation
correlations <- matrix(NA, nrow = nrow(X_CAIDE2), ncol = ncol(factors))
for (i in 1:ncol(factors)) {
  correlations[,i] <- apply(X_CAIDE2, 1, 
                            function(x){cor(x, 
                                            factors[,i], 
                                            method = "spearman")})
}
rownames(correlations) <- rownames(X_CAIDE2)
colnames(correlations) <- colnames(factors)
save(correlations, file = paste0(DataDir,"correlations_CAIDE2factors.RData"))

# Select top correlated features for each factor
probes <- list()
for (p in 1:7){
  probes[[p]] <- names(tail(sort(abs(correlations[,p])),1600))
}

# get exactly 10,000 probes
n = 1
finalProbes <- unique(unlist(probes))
while (length(finalProbes) > 10000){
  probes[[n]] <- probes[[n]][-1]
  finalProbes <- unique(unlist(probes))
  
  if (n < 7){
    n = n + 1
  } else {
    n = 1
  }
}

X_nonTest_Cor2 <- X_nonTest[finalProbes, ]
X_CAIDE2_Cor2 <- X_CAIDE2[finalProbes, ]
X_test_Cor2 <- X_test[finalProbes, ]

save(X_nonTest_Cor2, file = "CAIDE2_Cor/X_nonTest_Cor2.RData")
save(X_CAIDE1_Cor2, file = "CAIDE2_Cor/X_CAIDE2_Cor2.RData")
save(X_test_Cor2, file = "CAIDE2_Cor/X_test_Cor2.RData")

# Probes for CV
CVProbes <- list()
for (j in colnames(correlations)){
  CVProbes[[j]] <- names(tail(sort(abs(correlations[,j])),10000))
}
CVProbes <- unique(unlist(CVProbes))

X_CAIDE2_Cor2CV <- X_CAIDE2[CVProbes, ]
save(X_CAIDE2_Cor2CV, file = "CAIDE1_Cor/X_CAIDE2_Cor2CV.RData")


