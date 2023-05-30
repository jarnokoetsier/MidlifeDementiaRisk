# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

#BiocManager::install("sva")
library(sva)

# Load data
load("Data/X_nonTest.RData")
load("ADNI/X_ADNI.RData")

# combine data
dataMatrix <- cbind(X_nonTest, X_ADNI)

# convert to M values
dataMatrix_M <- log2(dataMatrix/(1 - dataMatrix))

# Create pheno data matrix
pheno <- data.frame(ID = c(colnames(X_nonTest), colnames(X_ADNI)),
                    batch = c(rep("EXTEND", ncol(X_nonTest)),
                              rep("ADNI", ncol(X_ADNI))))

# Model
modcombat <- model.matrix(~1, data=pheno)

# Batch correction
adj_DataMatrix <- ComBat(
  dataMatrix_M,
  batch,
  mod = modcombat,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = "EXTEND",
  BPPARAM = bpparam("SerialParam")
)