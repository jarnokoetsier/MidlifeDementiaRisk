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

# Make training and test data
load(paste0(DataDir,"CAIDE.Rdata"))
load(paste0(DataDir,"CAIDE2.Rdata"))
load(paste0(DataDir,"EPILIBRA.Rdata"))
load(paste0(DataDir,"metaData_ageFil.RData"))
load(paste0(DataDir,"methSet_allNorm_fil.RData"))
load(paste0(DataDir,"cellType.RData"))
load(paste0(DataDir,"TestTrain.RData"))

# Load machine learning functions
source("FUN_MachineLearning.R")

# Remove X and Y chromosomal genes
load("probe_annotation.RData")
probe_ann_fil <- probe_annotation[(probe_annotation$Chr != "chrX") &
                                    (probe_annotation$Chr != "chrY"), ]
methSet_allNorm_fil <- methSet_allNorm_fil[rownames(methSet_allNorm_fil) %in% probe_ann_fil$ID,]

# All
all_X <- as.matrix(t(methSet_allNorm_fil))
all_Y <- dat[,c("ID", "Basename", "Position", "Plate", "Age", "Sex")]
all_Y <- inner_join(all_Y, cellType, by = c("Basename" = "ID"))
all(rownames(all_X) == all_Y$Basename)

# Non-test
length(intersect(all_Y$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_nonTest <- t(all_X[setdiff(all_Y$Basename, TestTrain$Test),])
Y_nonTest <- all_Y[all_Y$Basename %in% setdiff(all_Y$Basename, TestTrain$Test),]
all(colnames(X_nonTest) == Y_nonTest$Basename)

# Save
save(X_nonTest, file = paste0(DataDir,"X_nonTest.RData"))
save(Y_nonTest, file = paste0(DataDir,"Y_nonTest.RData"))

# Test
X_test <- t(all_X[TestTrain$Test,])
Y_test <- all_Y[all_Y$Basename %in% TestTrain$Test,]
Y_test <- inner_join(Y_test, CAIDE[,c(1,17)], by = c("ID" = "ID"))
Y_test <- inner_join(Y_test, CAIDE2[,c(1,18)], by = c("ID" = "ID"))
Y_test <- inner_join(Y_test, EPILIBRA[,c(1,21)], by = c("ID" = "ID"))
X_test <- X_test[,Y_test$Basename]
all(colnames(X_test) == Y_test$Basename)

# Save
save(X_test, file = paste0(DataDir,"X_test.RData"))
save(Y_test, file = paste0(DataDir,"Y_test.RData"))

#*****************************************************************************#
# CAIDE1
#*****************************************************************************#

# CAIDE1: Training
length(intersect(CAIDE$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_CAIDE1 <- t(all_X[setdiff(CAIDE$Basename, TestTrain$Test),])
Y_CAIDE1 <- all_Y[all_Y$Basename %in% setdiff(CAIDE$Basename, TestTrain$Test),]
Y_CAIDE1 <- inner_join(Y_CAIDE1, CAIDE[,c(1,10:17)], by = c("ID" = "ID"))
all(colnames(X_CAIDE1) == Y_CAIDE1$Basename)

# Save
save(X_CAIDE1, file = paste0(DataDir,"X_CAIDE1.RData"))
save(Y_CAIDE1, file = paste0(DataDir,"Y_CAIDE1.RData"))

# Create cross-validation Index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(Y_CAIDE1),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)

# Load CV index
save(CVindex,file = paste0(DataDir,"CVindex_CAIDE1.RData"))

#*****************************************************************************#
# CAIDE2
#*****************************************************************************#

# CAIDE2: Training
length(intersect(CAIDE2$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_CAIDE2 <- t(all_X[setdiff(CAIDE2$Basename, TestTrain$Test),])
Y_CAIDE2 <- all_Y[all_Y$Basename %in% setdiff(CAIDE2$Basename, TestTrain$Test),]
Y_CAIDE2 <- inner_join(Y_CAIDE2, CAIDE2[,c(1,10:18)], by = c("ID" = "ID"))
all(colnames(X_CAIDE2) == Y_CAIDE2$Basename)

# Save
save(X_CAIDE2, file = paste0(DataDir,"X_CAIDE2.RData"))
save(Y_CAIDE2, file = paste0(DataDir,"Y_CAIDE2.RData"))

# Create cross-validation Index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(Y_CAIDE2),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)

# Load CV index
save(CVindex,file = paste0(DataDir,"CVindex_CAIDE2.RData"))


#*****************************************************************************#
# LIBRA
#*****************************************************************************#
# LIBRA: Training
length(intersect(EPILIBRA$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_LIBRA <- t(all_X[setdiff(EPILIBRA$Basename, TestTrain$Test),])
Y_LIBRA <- all_Y[all_Y$Basename %in% setdiff(EPILIBRA$Basename, TestTrain$Test),]
Y_LIBRA <- inner_join(Y_LIBRA, EPILIBRA[,c(1,10:21)], by = c("ID" = "ID"))
all(colnames(X_LIBRA) == Y_LIBRA$Basename)

# Save
save(X_LIBRA, file = paste0(DataDir,"X_LIBRA.RData"))
save(Y_LIBRA, file = paste0(DataDir,"Y_LIBRA.RData"))

# Create cross-validation Index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(Y_LIBRA),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)

# Load CV index
save(CVindex,file = paste0(DataDir,"CVindex_LIBRA.RData"))
