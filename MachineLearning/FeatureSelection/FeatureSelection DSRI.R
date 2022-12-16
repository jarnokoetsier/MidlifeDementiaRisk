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
source("FUN_MachineLearning.R")

#*****************************************************************************#
# Make all data combinations
#*****************************************************************************#
# Make training and test data
load("CAIDE.Rdata")
load("CAIDE2.Rdata")
load("EPILIBRA.Rdata")
load("metaData_ageFil.RData")
load("methSet_allNorm_fil.RData")
load("cellType.RData")
load("gt_results.RData")
load("TestTrain.RData")

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

# CAIDE1
length(intersect(CAIDE$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_CAIDE1 <- t(all_X[setdiff(CAIDE$Basename, TestTrain$Test),])
Y_CAIDE1 <- all_Y[all_Y$Basename %in% setdiff(CAIDE$Basename, TestTrain$Test),]
Y_CAIDE1 <- inner_join(Y_CAIDE1, CAIDE[,c(1,10:17)], by = c("ID" = "ID"))
all(colnames(X_CAIDE1) == Y_CAIDE1$Basename)

# CAIDE2
length(intersect(CAIDE2$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_CAIDE2 <- t(all_X[setdiff(CAIDE2$Basename, TestTrain$Test),])
Y_CAIDE2 <- all_Y[all_Y$Basename %in% setdiff(CAIDE2$Basename, TestTrain$Test),]
all(colnames(X_CAIDE2) == Y_CAIDE2$Basename)

# LIBRA
length(intersect(EPILIBRA$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_LIBRA <- t(all_X[setdiff(EPILIBRA$Basename, TestTrain$Test),])
Y_LIBRA <- all_Y[all_Y$Basename %in% setdiff(EPILIBRA$Basename, TestTrain$Test),]
Y_LIBRA <- inner_join(Y_LIBRA, EPILIBRA[,c(1,10:21)], by = c("ID" = "ID"))
all(colnames(X_LIBRA) == Y_LIBRA$Basename)

# Non-test
length(intersect(all_Y$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_nonTest <- t(all_X[setdiff(all_Y$Basename, TestTrain$Test),])
Y_nonTest <- all_Y[all_Y$Basename %in% setdiff(all_Y$Basename, TestTrain$Test),]
all(colnames(X_nonTest) == Y_nonTest$Basename)

# Test
X_test <- t(all_X[TestTrain$Test,])
Y_test <- all_Y[all_Y$Basename %in% TestTrain$Test,]
Y_test <- inner_join(Y_test, EPILIBRA[,c(1,10:21)], by = c("ID" = "ID"))
Y_test <- inner_join(Y_test, CAIDE[,c(1,10:17)], by = c("ID" = "ID"))

X_test <- X_test[,Y_test$Basename]
all(colnames(X_test) == Y_test$Basename)

#*****************************************************************************#
# No selection
#*****************************************************************************#

# Get beta- and M-values
Mvalues <- log2(X_nonTest/(1 - X_nonTest))

# Format data
plot_all <- data.frame(
  ID = rownames(X_nonTest),
  meanBeta = rowMeans(X_nonTest),
  sdBeta = apply(X_nonTest, 1, sd),
  meanM  = rowMeans(Mvalues),
  sdM = apply(Mvalues, 1, sd)
)

save(plot_all, file = "plot_all.RData")
#*****************************************************************************#
# Variance selection
#*****************************************************************************#

cpg_var <- apply(X_nonTest, 1, var)
cpg_selected_var <- names(tail(sort(cpg_var), 10000))
rm(cpg_var)
#all(rownames(X_nonTest) %in% probe_annotation$ID)

X_nonTest_var <- X_nonTest[cpg_selected_var, ]
X_CAIDE1_var <- X_CAIDE1[cpg_selected_var, ]
X_LIBRA_var <- X_LIBRA[cpg_selected_var, ]
X_test_var <- X_test[cpg_selected_var, ]

save(X_nonTest_var, file = "X_nonTest_var.RData")
save(X_CAIDE1_var, file = "X_CAIDE1_var.RData")
save(X_LIBRA_var, file = "X_LIBRA_var.RData")
save(X_test_var, file = "X_test_var.RData")

rm(X_nonTest_var)
rm(X_CAIDE1_var)
rm(X_LIBRA_var)
rm(X_test_var)

save(Y_nonTest, file = "Y_nonTest.RData")
save(Y_CAIDE1, file = "Y_CAIDE1.RData")
save(Y_LIBRA, file = "Y_LIBRA.RData")
save(Y_test, file = "Y_test.RData")

#*****************************************************************************#
# Variance selection (based on M-value)
#*****************************************************************************#

X_nonTest_M <- log2(X_nonTest/(1-X_nonTest))
cpg_var <- apply(X_nonTest_M, 1, var)
cpg_selected_var <- names(tail(sort(cpg_var), 10000))
rm(cpg_var)
#all(rownames(X_nonTest) %in% probe_annotation$ID)

X_nonTest_varM <- X_nonTest[cpg_selected_var, ]
X_CAIDE1_varM <- X_CAIDE1[cpg_selected_var, ]
X_LIBRA_varM <- X_LIBRA[cpg_selected_var, ]
X_test_varM <- X_test[cpg_selected_var, ]

save(X_nonTest_varM, file = "X_nonTest_varM.RData")
save(X_CAIDE1_varM, file = "X_CAIDE1_varM.RData")
save(X_LIBRA_varM, file = "X_LIBRA_varM.RData")
save(X_test_varM, file = "X_test_varM.RData")


#*****************************************************************************#
# S-score selection
#*****************************************************************************#

calculate_S <- function(x){
  S = abs(mean(x)- 0.5)/var(x)
}

cpg_S <- apply(X_nonTest, 1, calculate_S)
cpg_selected_S <- names(tail(sort(cpg_S),10000))
rm(cpg_S)

X_nonTest_S <- X_nonTest[cpg_selected_S, ]
X_CAIDE1_S <- X_CAIDE1[cpg_selected_S, ]
X_LIBRA_S <- X_LIBRA[cpg_selected_S, ]
X_test_S <- X_test[cpg_selected_S, ]

save(X_nonTest_S, file = "X_nonTest_S.RData")
save(X_CAIDE1_S, file = "X_CAIDE1_S.RData")
save(X_LIBRA_S, file = "X_LIBRA_S.RData")
save(X_test_S, file = "X_test_S.RData")

rm(X_nonTest_S)
rm(X_CAIDE1_S)
rm(X_LIBRA_S)
rm(X_test_S)


#*****************************************************************************#
# KS-like selection
#*****************************************************************************#


# Convert to M-values
X_nonTest_M <- log2(X_nonTest/(1-X_nonTest))

# Get three bins
t1 <- quantile(plot_all$sdM,0.33)
t2 <- quantile(plot_all$sdM,0.67)

# Get groups: for each of these groups were are going the select the
# most representative features
probe_annotation$Relation_to_Island[(probe_annotation$Relation_to_Island != "Island") &
                                      (probe_annotation$Relation_to_Island != "OpenSea")] <- "ShelforShore"
Group <- data.frame(CpG = probe_annotation$ID,
                    Group = paste(probe_annotation$Class,
                                  probe_annotation$Relation_to_Island, 
                                  probe_annotation$Type, sep = "_"))

Group <- inner_join(Group, plot_all, by = c("CpG" = "ID"))
Group$Bins <- rep("Intermediate", nrow(Group))
Group$Bins[Group$sdM < t1] <- "Low"
Group$Bins[Group$sdM > t2] <- "High"
Group$Group1 <- paste0(Group$Group, "_", Group$Bins)
save(Group, file = "Group.RData")

# Number of features to select from each group
n <- 20000
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

# Make clusters
nCores = 54
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Record start time
t_start <- Sys.time()

# For each Group we are going to perform feature selection:
output <- foreach(i =  1:length(uniqueGroups), .inorder = FALSE) %dopar% {
  
  # Filter data for CpG in group
  dataMatrix_copy_fil <- dataMatrix_copy[rownames(dataMatrix_copy) %in% Group$CpG[Group$Group1 == uniqueGroups[i]],]   
  
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

save(probes, file = "probes.RData")





randomFeatures <- sample(rownames(X_nonTest), 80000)
save(randomFeatures, file = "randomFeatures.RData")

# Get features
load("randomFeatures.RData")
load("plot_all.RData")
plot_all <- plot_all[plot_all$ID %in% randomFeatures,]

# Convert to M-values
X_nonTest_M <- log2(X_nonTest/(1-X_nonTest))
#X_nonTest_M <- X_nonTest_M[randomFeatures,]

# Get three bins
t1 <- quantile(plot_all$sdM,0.33)
t2 <- quantile(plot_all$sdM,0.67)

# Get groups: for each of these groups were are going the select the
# most representative features
probe_annotation$Relation_to_Island[(probe_annotation$Relation_to_Island != "Island") &
                                      (probe_annotation$Relation_to_Island != "OpenSea")] <- "ShelforShore"
Group <- data.frame(CpG = probe_annotation$ID,
                    Group = paste(probe_annotation$Class,
                                  probe_annotation$Relation_to_Island, 
                                  probe_annotation$Type, sep = "_"))

Group <- inner_join(Group, plot_all, by = c("CpG" = "ID"))
Group$Bins <- rep("Intermediate", nrow(Group))
Group$Bins[Group$sdM < t1] <- "Low"
Group$Bins[Group$sdM > t2] <- "High"
Group$Group1 <- paste0(Group$Group, "_", Group$Bins)
save(Group, file = "Group.RData")

# Number of features to select from each group
n <- 20000
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


# Make clusters
nCores = 10
cl <- makeCluster(nCores)
registerDoParallel(cl)

# For each Group we are going to perform feature selection:
output <- foreach(i =  1:length(uniqueGroups), .combine = c, .inorder = FALSE) %do% {
  
  # Filter data for CpG in group
  dataMatrix_copy_fil <- dataMatrix_copy[rownames(dataMatrix_copy) %in% Group$CpG[Group$Group1 == uniqueGroups[i]],]   
  
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

save(output, file = "output_all.RData")

probes54 <- probes
save(probes54, file = "probes54.RData")

# From these 20,000 probes, we can now select 10,000 probes
load("output.RData")
test <- unlist(output)

# Select seed probe
cpg_var <- apply(X_nonTest_M[test,], 1, var)
seedProbe <- names(cpg_var)[which.max(cpg_var)]

# Perform KS-like feature selection
probes <- selectionKS(dataMatrix =  X_nonTest_M[test,],
                      nFeatures = 1000,
                      seedProbe = seedProbe)

save(probes, file = "probes.RData")
