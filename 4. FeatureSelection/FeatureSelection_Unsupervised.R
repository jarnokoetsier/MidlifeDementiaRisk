
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

# Set direcotries
DataDir <- "Data/"

# Load data:

# Data matrices:
load(paste0(DataDir,"X_CAIDE1.RData"))
load(paste0(DataDir,"X_CAIDE2.RData"))
load(paste0(DataDir,"X_LIBRA.RData"))
load(paste0(DataDir,"X_nonTest.RData"))
load(paste0(DataDir,"X_test_Cor.RData"))

# Target variables:
load(paste0(DataDir,"Y_CAIDE1.RData"))
load(paste0(DataDir,"Y_CAIDE2.RData"))
load(paste0(DataDir,"Y_LIBRA.RData"))
load(paste0(DataDir,"Y_nonTest.RData"))
load(paste0(DataDir,"Y_test.RData"))

# Test if samples are in correct order
all(colnames(X_train) == Y_CAIDE1$Basename)


###############################################################################

# Variance-based feature selection

###############################################################################

#*****************************************************************************#
# Beta-values
#*****************************************************************************#

# Select features
cpg_var <- apply(X_nonTest, 1, var)
cpg_selected_var <- names(tail(sort(cpg_var), 10000))
rm(cpg_var)

# Extract matrix
X_nonTest_var <- X_nonTest[cpg_selected_var, ]
X_CAIDE1_var <- X_CAIDE1[cpg_selected_var, ]
X_CAIDE2_var <- X_CAIDE2[cpg_selected_var, ]
X_LIBRA_var <- X_LIBRA[cpg_selected_var, ]
X_test_var <- X_test[cpg_selected_var, ]

# Save matrix
save(X_nonTest_var, file = "X_nonTest_var.RData")
save(X_CAIDE1_var, file = "X_CAIDE1_var.RData")
save(X_CAIDE2_var, file = "X_CAIDE2_var.RData")
save(X_LIBRA_var, file = "X_LIBRA_var.RData")
save(X_test_var, file = "X_test_var.RData")

# Clear workspace
rm(X_nonTest_var)
rm(X_CAIDE1_var)
rm(X_CAIDE2_var)
rm(X_LIBRA_var)
rm(X_test_var)

#*****************************************************************************#
# M-values
#*****************************************************************************#

# Convert Beta to M-values
X_nonTest_M <- log2(X_nonTest/(1-X_nonTest))

# Select features
cpg_var <- apply(X_nonTest_M, 1, var)
cpg_selected_var <- names(tail(sort(cpg_var), 10000))
rm(cpg_var)

# Extract matrix
X_nonTest_varM <- X_nonTest[cpg_selected_var, ]
X_CAIDE1_varM <- X_CAIDE1[cpg_selected_var, ]
X_CAIDE2_varM <- X_CAIDE2[cpg_selected_var, ]
X_LIBRA_varM <- X_LIBRA[cpg_selected_var, ]
X_test_varM <- X_test[cpg_selected_var, ]

# Save matrix
save(X_nonTest_varM, file = "X_nonTest_varM.RData")
save(X_CAIDE1_varM, file = "X_CAIDE1_varM.RData")
save(X_CAIDE2_varM, file = "X_CAIDE2_varM.RData")
save(X_LIBRA_varM, file = "X_LIBRA_varM.RData")
save(X_test_varM, file = "X_test_varM.RData")

# Clear workspace
rm(X_nonTest_varM)
rm(X_CAIDE1_varM)
rm(X_CAIDE2_varM)
rm(X_LIBRA_varM)
rm(X_test_varM)


#*****************************************************************************#
# M-values (corrected)
#*****************************************************************************#

# Convert beta to M
X_nonTest_M <- log2(X_nonTest/(1-X_nonTest))

# Get Features
features <- rownames(X_nonTest_M)

# Make formula
formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                 paste(colnames(Y_nonTest[,7:12]), collapse = " + "))

# Combine with cell type composition
dataMatrix_scaled <- cbind(dataMatrix_scaled[,features],Y_nonTest[,7:12])

# Make linear model
model <- lm(as.formula(formula), data = as.data.frame(dataMatrix_scaled))

# Get residual values
residualValues <- residuals(model)
save(residualValues, file = "residualValues.RData")

res <- NULL
for (i in 1:8){
  load(paste0("residualValues",i,".RData"))
  res <- cbind(res, residualValues)
}

# Save residual values
save(res, file = "res_all.RData")

# Select features
cpg_var <- apply(res, 2, var)
cpg_selected_varCor <- names(tail(sort(cpg_var), 10000))
rm(cpg_var)

# Extract matrix
X_nonTest_varMCor <- X_nonTest[cpg_selected_varCor, ]
X_CAIDE1_varMCor <- X_CAIDE1[cpg_selected_varCor, ]
X_CAIDE2_varMCor <- X_CAIDE2[cpg_selected_varCor, ]
X_LIBRA_varMCor <- X_LIBRA[cpg_selected_varCor, ]
X_test_varMCor <- X_test[cpg_selected_varCor, ]

# Save Matrix
save(X_nonTest_varMCor, file = "X_nonTest_varMCor.RData")
save(X_CAIDE1_varMCor, file = "X_CAIDE1_varMCor.RData")
save(X_CAIDE2_varMCor, file = "X_CAIDE2_varMCor.RData")
save(X_LIBRA_varMCor, file = "X_LIBRA_varMCor.RData")
save(X_test_varMCor, file = "X_test_varMCor.RData")

# Clear workspace
rm(X_nonTest_varMCor)
rm(X_CAIDE1_varMCor)
rm(X_CAIDE2_varMCor)
rm(X_LIBRA_varMCor)
rm(X_test_varMCor)

#*****************************************************************************#
# Beta-values (corrected)
#*****************************************************************************#

# Convert residual M values to residual Beta values
resB <- (2^res)/((2^res) + 1)

# Select features
cpg_var <- apply(resB, 2, var)
cpg_selected_varCor <- names(tail(sort(cpg_var), 10000))
rm(cpg_var)

# Extract matrix
X_nonTest_varCor <- X_nonTest[cpg_selected_varCor, ]
X_CAIDE1_varCor <- X_CAIDE1[cpg_selected_varCor, ]
X_CAIDE2_varCor <- X_CAIDE2[cpg_selected_varCor, ]
X_LIBRA_varCor <- X_LIBRA[cpg_selected_varCor, ]
X_test_varCor <- X_test[cpg_selected_varCor, ]

# Save matrix
save(X_nonTest_varCor, file = "X_nonTest_varCor.RData")
save(X_CAIDE1_varCor, file = "X_CAIDE1_varCor.RData")
save(X_CAIDE2_varCor, file = "X_CAIDE2_varCor.RData")
save(X_LIBRA_varCor, file = "X_LIBRA_varCor.RData")
save(X_test_varCor, file = "X_test_varCor.RData")

# Clear workspace
rm(X_nonTest_varCor)
rm(X_CAIDE1_varCor)
rm(X_CAIDE2_varCor)
rm(X_LIBRA_varCor)
rm(X_test_varCor)



###############################################################################

# S-score-based feature selection

###############################################################################

# Formula of S-score
calculate_S <- function(x){
  S = abs(mean(x)- 0.5)/var(x)
}

# select features
cpg_S <- apply(X_nonTest, 1, calculate_S)
cpg_selected_S <- names(tail(sort(cpg_S),10000))
rm(cpg_S)

# Extract matrix
X_nonTest_S <- X_nonTest[cpg_selected_S, ]
X_CAIDE1_S <- X_CAIDE1[cpg_selected_S, ]
X_CAIDE2_S <- X_CAIDE2[cpg_selected_S, ]
X_LIBRA_S <- X_LIBRA[cpg_selected_S, ]
X_test_S <- X_test[cpg_selected_S, ]

# save matrix
save(X_nonTest_S, file = "X_nonTest_S.RData")
save(X_CAIDE1_S, file = "X_CAIDE1_S.RData")
save(X_CAIDE2_S, file = "X_CAIDE2_S.RData")
save(X_LIBRA_S, file = "X_LIBRA_S.RData")
save(X_test_S, file = "X_test_S.RData")

# Clear workspace
rm(X_nonTest_S)
rm(X_CAIDE1_S)
rm(X_CAIDE2_S)
rm(X_LIBRA_S)
rm(X_test_S)


###############################################################################

# PCA-based feature selection

###############################################################################

# Perform PCA
pcaList <-  prcomp(t(X_nonTest),        
                   retx = TRUE,
                   center =TRUE,
                   scale = TRUE)

# Get loadings
loadings <- pcaList$rotation
X_nonTest_PC <- pcaList$x

# Get PC features
X_CAIDE1_scaled <- (X_CAIDE1 - rowMeans(X_CAIDE1))/(apply(X_CAIDE1, 1, sd))
X_CAIDE1_PC <- t(X_CAIDE1_scaled) %*% loadings

X_CAIDE2_scaled <- (X_CAIDE2 - rowMeans(X_CAIDE2))/(apply(X_CAIDE2, 1, sd))
X_CAIDE2_PC <- t(X_CAIDE2_scaled) %*% loadings

X_LIBRA_scaled <- (X_LIBRA - rowMeans(X_LIBRA))/(apply(X_LIBRA, 1, sd))
X_LIBRA_PC <- t(X_LIBRA_scaled) %*% loadings

X_test_scaled <- (X_test - rowMeans(X_test))/(apply(X_test, 1, sd))
X_test_PC <- t(X_test_scaled) %*% loadings

# Save matrix
save(X_nonTest_PC, file = "X_nonTest_PC.RData")
save(X_CAIDE1_PC, file = "X_CAIDE1_PC.RData")
save(X_CAIDE2_PC, file = "X_CAIDE2_PC.RData")
save(X_LIBRA_PC, file = "X_LIBRA_PC.RData")
save(X_test_PC, file = "X_test_PC.RData")

# Clear workspace
rm(X_nonTest_PC)
rm(X_CAIDE1_PC)
rm(X_CAIDE2_PC)
rm(X_LIBRA_PC)
rm(X_test_PC)

###############################################################################

# KS-like feature selection

###############################################################################

# Convert to M-values
X_nonTest_M <- log2(X_nonTest/(1-X_nonTest))


# Get groups: 
# for each of these groups were are going the select the most representative features:

# Get three bins based on the standard deviation
t1 <- quantile(plot_all$sdM,0.33)
t2 <- quantile(plot_all$sdM,0.67)

# Island vs open sea vs shelf or shore
probe_annotation$Relation_to_Island[(probe_annotation$Relation_to_Island != "Island") &
                                      (probe_annotation$Relation_to_Island != "OpenSea")] <- "ShelforShore"
# Make groups
Group <- data.frame(CpG = probe_annotation$ID,
                    Group = paste(probe_annotation$Class,
                                  probe_annotation$Relation_to_Island, 
                                  probe_annotation$Type, sep = "_"))

Group <- inner_join(Group, plot_all, by = c("CpG" = "ID"))
Group$Bins <- rep("Intermediate", nrow(Group))
Group$Bins[Group$sdM < t1] <- "Low"
Group$Bins[Group$sdM > t2] <- "High"
Group$Group1 <- paste0(Group$Group, "_", Group$Bins)

# Save probe groups
save(Group, file = "Group.RData")

# Calculate the number of features to select from each group
n <- 20000 # Total number of selected features
uniqueGroups <- unique(Group$Group1)
nFeatures <- rep(NA, length(uniqueGroups))
for (i in 1:length(uniqueGroups)){
  nFeatures[i] <- round(n*(table(Group$Group1)[uniqueGroups[i]]/sum(table(Group$Group1))))
}
nFeatures[which.max(nFeatures)] <- nFeatures[which.max(nFeatures)]-(n - sum(nFeatures))
names(nFeatures) <- uniqueGroups
nFeatures

# Make copy of data
dataMatrix_copy <- X_nonTest_M

# remove object
rm(all_X)
rm(all_Y)
rm(cellType)
rm(dat)
rm(methSet_allNorm_fil)
rm(plot_all)
rm(probe_ann_fil)
rm(probe_annotation)
rm(X_nonTest)
rm(X_nonTest_M)
rm(Y_nonTest)
rm(i)
rm(n)
rm(t1)
rm(t2)
gc()

nFeatures1 <- nFeatures


output <- NULL
for (i in 1:length(uniqueGroups)){
  
  # Filter data for CpG in group
  dataMatrix_copy_fil <- dataMatrix_copy[rownames(dataMatrix_copy) %in% Group$CpG[Group$Group1 == uniqueGroups[i]],]   
  
  # Seed probe (most variable probe)
  cpg_var <- apply(dataMatrix_copy_fil, 1, var)
  seedProbe <- names(cpg_var)[which.max(cpg_var)]
  
  # Make copy of data
  dataMatrix = dataMatrix_copy_fil
  nFeatures = nFeatures1[i]
  
  # Make matrix to save correlations
  cor_matrix <- matrix(NA, nrow(dataMatrix), nFeatures-1)
  rownames(cor_matrix) <- rownames(dataMatrix)
  
  # Make vector to save feature set
  newProbe <- rep(NA, nFeatures)
  newProbe[1] <- seedProbe
  
  # Start selecting features
  for (j in 1:nFeatures){
    
    # Calculate correlations between seed probe and all other probes
    cor_matrix[,j] <- apply(dataMatrix,1,function(x){cor(x,dataMatrix[newProbe[j],],
                                                         method = "spearman")})
    
    # Add most uncorrelated probe to feature set
    temp <- apply(as.data.frame(abs(cor_matrix[,1:j])), 1, max)
    newProbe[j+1] <- rownames(cor_matrix)[which.min(temp)]
    
    # Remove all highly correlated probes: keep probes with cor < 0.5
    cor_matrix <- cor_matrix[abs(cor_matrix[,j]) < 0.8,]
    dataMatrix <- dataMatrix[rownames(cor_matrix),]
  }
  
  output <- c(output,newProbe)
  save(output, file = "probesOutput.RData")
}

# Perform again on the 20,000 selected probes
dataMatrix_copy <- X_nonTest_M[output, ]

# Seed probe (most variable probe)
cpg_var <- apply(dataMatrix_copy_fil, 1, var)
seedProbe <- names(cpg_var)[which.max(cpg_var)]

# Make copy of data
dataMatrix = dataMatrix_copy_fil
nFeatures = nFeatures1[i]

# Make matrix to save correlations
cor_matrix <- matrix(NA, nrow(dataMatrix), nFeatures-1)
rownames(cor_matrix) <- rownames(dataMatrix)

# Make vector to save feature set
newProbe <- rep(NA, nFeatures)
newProbe[1] <- seedProbe

# Start selecting features
for (j in 1:10000){
  
  # Calculate correlations between seed probe and all other probes
  cor_matrix[,j] <- apply(dataMatrix,1,function(x){cor(x,dataMatrix[newProbe[j],],
                                                       method = "spearman")})
  
  # Add most uncorrelated probe to feature set
  temp <- apply(as.data.frame(abs(cor_matrix[,1:j])), 1, max)
  newProbe[j+1] <- rownames(cor_matrix)[which.min(temp)]
  
  # Remove all highly correlated probes: keep probes with cor < 0.5
  cor_matrix <- cor_matrix[abs(cor_matrix[,j]) < 0.8,]
  dataMatrix <- dataMatrix[rownames(cor_matrix),]
}

probes<- newProbe
save(probes, file = "probesOutput_final.RData")


# Select features from data
X_nonTest_KS <- X_nonTest_KS[probes, ]
X_CAIDE1_KS <- X_CAIDE1[probes, ]
X_CAIDE2_KS <- X_CAIDE2[probes, ]
X_LIBRA_KS <- X_LIBRA[probes, ]
X_test_KS <- X_test[probes, ]

# Save data
save(X_nonTest_KS, file = "X_nonTest_KS.RData")
save(X_CAIDE1_KS, file = "X_CAIDE1_KS.RData")
save(X_CAIDE2_KS, file = "X_CAIDE2_KS.RData")
save(X_LIBRA_KS, file = "X_LIBRA_KS.RData")
save(X_test_KS, file = "X_test_KS.RData")


