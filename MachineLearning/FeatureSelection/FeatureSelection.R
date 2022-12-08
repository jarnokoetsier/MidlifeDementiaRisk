# Load packages
library(doParallel)
library(foreach)
library(ggvenn)
library(glmnet)
library(caret)

# Load machine learning functions
source("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA/MachineLearning/FUN_MachineLearning.R")

# Set working directory
setwd("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA")

# Load normalized training data
load("E:/Thesis/MLData/mydat1.RData")
load("E:/Thesis/MLData/pheno.RData")
dataMatrix <- unique(mydat1)

###############################################################################

# Most variable features (Beta-values)

###############################################################################

# Calculate variance
cpg_var <- apply(dataMatrix, 1, var)

# Get 10,000 probes with highest variance
cpg_selected_var <- names(tail(sort(cpg_var), 10000))

# Select the probes in the data matrix
dataMatrix_var <- dataMatrix[cpg_selected_var, ]

# Save output
save(dataMatrix_var, file = "dataMatrix_var.RData")

###############################################################################

# Most variable features (M-values)

###############################################################################

# Calculate variance
cpg_varM <- apply(dataMatrix_M, 1, var)

# Get 10,000 probes with highest variance
cpg_selected_varM <- names(tail(sort(cpg_varM), 10000))

# Select the probes in the data matrix
dataMatrix_varM <- dataMatrix[cpg_selected_varM, ]

# Save output
save(dataMatrix_varM, file = "dataMatrix_var.RData")


###############################################################################

# Features with highest S-score

###############################################################################

# Function of S-score
calculate_S <- function(x){
  S = abs(mean(x)- 0.5)/var(x)
}

# Calculate S-score
cpg_S <- apply(dataMatrix, 1, calculate_S)

# Get 10,000 probes with highest S-score
cpg_selected_S <- names(tail(sort(cpg_S),10000))


###############################################################################

# Kennard-Stone algorithm

###############################################################################

# Convert beta-values to M-values
dataMatrix_fil <- dataMatrix[rowSums((dataMatrix > 0) & (dataMatrix < 1)) == ncol(dataMatrix), ]
dataMatrix_M <- log2(dataMatrix_fil/(1 - dataMatrix_fil))

# Get groups: for each of these groups were are going the select the
# most representative features

# Load annotation data
load("E:/Thesis/MLData/probe_annotation.RData")

# Location with respect to CpG islands
probe_annotation$Relation_to_Island[(probe_annotation$Relation_to_Island != "Island") &
                                      (probe_annotation$Relation_to_Island != "OpenSea")] <- "ShelforShore"

# 18 groups based on location with respect to genes, location with respect to
# CpG islands, and probe type
Group <- data.frame(CpG = probe_annotation$ID,
                    Group = paste(probe_annotation$Class,
                                  probe_annotation$Relation_to_Island, 
                                  probe_annotation$Type, sep = "_"))

Group <- Group[Group$CpG %in% rownames(dataMatrix_M),]

# Calculate the number of features to select from each group
# The larger the group, the more features to select
n <- 5000
uniqueGroups <- unique(Group$Group)
nFeatures <- rep(NA, length(uniqueGroups))
for (i in 1:length(uniqueGroups)){
  nFeatures[i] <- round(n*(table(Group$Group)[uniqueGroups[i]]/sum(table(Group$Group))))
}
nFeatures[which.max(nFeatures)] <- nFeatures[which.max(nFeatures)]-(n - sum(nFeatures))
names(nFeatures) <- uniqueGroups

# Make copy of data
dataMatrix_copy <- dataMatrix_M

# Make clusters
nCores = 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Record start time
t_start <- Sys.time()

# For each Group we are going to perform feature selection:
output <- foreach(i =  1:length(uniqueGroups), .inorder = FALSE) %dopar% {
  
  # Filter data for CpG in group
  dataMatrix_copy_fil <- dataMatrix_copy[rownames(dataMatrix_copy) %in% Group$CpG[Group$Group == uniqueGroups[i]],]   
  
  # Seed probe (most variable probe)
  cpg_var <- apply(dataMatrix_copy_fil, 1, var)
  seedProbe <- names(cpg_var)[which.max(cpg_var)]
  
  # perform feature selection
  probes <- selectionKS(dataMatrix = dataMatrix_copy_fil,
              nFeatures = nFeatures[i],
              seedProbe = seedProbe)
  
  return(probes)
}
#Stop clusters
stopCluster(cl)

# Record end time
t_end <- Sys.time()

# Give run time
t_end-t_start


###############################################################################

# Unsupervised LASSO

###############################################################################

# Convert beta-values to M-values
dataMatrix_fil <- dataMatrix[rowSums((dataMatrix > 0) & (dataMatrix < 1)) == ncol(dataMatrix), ]
dataMatrix_M <- log2(dataMatrix_fil/(1 - dataMatrix_fil))

# Scale the data
dataMatrix_scaled <- (dataMatrix_M - rowMeans(dataMatrix_M))/(apply(dataMatrix_M, 1, sd))

# Perform permutation
switch <- function(x) {sample(x, length(x), replace = FALSE)}
dataMatrix_perm <- t(apply(dataMatrix_scaled,1,switch))
dataMatrix_perm <- apply(dataMatrix_perm,2,switch)
dataMatrix_all <- cbind(dataMatrix_scaled, dataMatrix_perm)


# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           number = 5, 
                           repeats = 5, 
                           search = "grid", 
                           savePredictions = FALSE)

# Set grid for lambda (2.5 for CAIDE)
lambdaCV <- exp(seq(log(0.001),log(1),length.out = 50))

# Set grid for alpha
alphaCV <- 1

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Use MSE as performance metric
performance_metric = "Accuracy"
MLmethod = "glmnet"

# Register cores for parallel computing
#detectCores()
#nCores <- 3
#cl <- makeCluster(nCores)
#registerDoParallel(cl)

# Actual training
set.seed(123)
fit <- train(x = t(dataMatrix_all), 
             y = factor(c(rep("Normal", ncol(dataMatrix_scaled)), 
                          rep("Permuted", ncol(dataMatrix_perm)))), 
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)

# Stop clusters
#stopCluster(cl)

test <- as.matrix(coef(fit$finalModel))
coeffs <- test[-1, which.min(abs(colSums(test != 0) - 10000))]
probes <- coeffs[coeffs != 0]

# Get results
trainResults <- fit$results

# Get optimal lambda and alpha
optAlpha <- fit$bestTune$alpha
optLambda <- fit$bestTune$lambda

finalModel <- glmnet(x = t(dataMatrix_all), 
                      y = factor(c(rep("Normal", ncol(dataMatrix_scaled)), 
                                   rep("Permuted", ncol(dataMatrix_perm)))), 
                      family = "binomial",
                      alpha = optAlpha, 
                      lambda = optLambda,
                      standardize = TRUE)

test <- as.matrix(coef(finalModel))




