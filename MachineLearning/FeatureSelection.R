library(doParallel)
library(foreach)

# Set working directory
setwd("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA")

# Load normalized training data
load("E:/Thesis/MLData/mydat1.RData")
load("E:/Thesis/MLData/pheno.RData")

#*****************************************************************************#
# Select most variable features
#*****************************************************************************#

cpg_var <- apply(mydat1, 1, var)
cpg_selected_var <- names(tail(sort(cpg_var), 5000))


#*****************************************************************************#
# Select most variable features
#*****************************************************************************#

calculate_S <- function(x){
  S = abs(mean(x)- 0.5)/var(x)
}

cpg_S <- apply(mydat1, 1, calculate_S)
cpg_selected_S <- names(tail(sort(cpg_S),5000))

length(intersect(cpg_selected_var, cpg_selected_S))

#*****************************************************************************#
# Kennard-Stone-like feature selection
#*****************************************************************************#
dataMatrix <- mydat1

# Convert to M-values
dataMatrix_fil <- mydat1[rowSums((dataMatrix > 0) & (dataMatrix < 1)) == ncol(dataMatrix), ]
dataMatrix_M <- log2(dataMatrix_fil/(1 + dataMatrix_fil))

# Number of features to select
nFeatures <- c(200,200,100)

# Groups
load("E:/Thesis/MLData/probe_annotation.RData")
Group <- data.frame(CpG = probe_annotation$ID,
                    Group = probe_annotation$Class)

# Make copy of data
dataMatrix_copy <- dataMatrix_M

# Make list for selected probes
selectedProbes <- list()

# Make clusters
nCores = 3
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Record start time
t_start <- Sys.time()

uniqueGroups <- unique(Group[,2])



# For each parameter combination....
output <- foreach(i =  1:length(uniqueGroups), .inorder = FALSE) %dopar% {
  
  # Filter data for Group
  dataMatrix_copy_fil <- dataMatrix_copy[rownames(dataMatrix_copy) %in% Group$CpG[Group$Group == uniqueGroups[i]],]   
  
  # Seed probe
  cpg_var <- apply(dataMatrix_copy_fil, 1, var)
  seedProbe <- names(cpg_var)[which.max(cpg_var)]
  
  # Make matrix to save correlations
  cor_matrix <- matrix(NA, nrow(dataMatrix_copy_fil), nFeatures[i]-1)
  rownames(cor_matrix) <- rownames(dataMatrix_copy_fil)
  
  # Make vector to save feature set
  newProbe <- rep(NA, nFeatures[i])
  newProbe[1] <- seedProbe
  
  # Start selecting features
  for (j in 1:(nFeatures[i] - 1)){
    # Calculate correlations between seed probe and all other probes
    cor_matrix[,j] <- apply(dataMatrix_copy_fil,1,function(x){cor(x,dataMatrix_copy_fil[newProbe[j],])})
    
    # Add most uncorrelated probe to feature set
    newProbe[j+1] <- rownames(cor_matrix)[which.min(abs(cor_matrix[,j]))]
    
    # Remove all highly correlated probes: keep probes with cor < 0.9
    cor_matrix <- cor_matrix[abs(cor_matrix[,j]) < 0.5,]
    dataMatrix_copy_fil <- dataMatrix_copy_fil[rownames(cor_matrix),]
  }
  return(newProbe)
}
#Stop clusters
stopCluster(cl)

# Record end time
t_end <- Sys.time()

# Give run time
t_end-t_start



for (i in 1:length(uniqueGroups)) {
 
  # Filter data for Group
  dataMatrix_copy_fil <- dataMatrix_copy[rownames(dataMatrix_copy) %in% Group$CpG[Group$Group == uniqueGroups[i]],]   
  
  # Seed probe
  cpg_var <- apply(dataMatrix_copy_fil, 1, var)
  seedProbe <- names(cpg_var)[which.max(cpg_var)]
  
  # Make matrix to save correlations
  cor_matrix <- matrix(NA, nrow(dataMatrix_copy_fil), nFeatures[i]-1)
  rownames(cor_matrix) <- rownames(dataMatrix_copy_fil)
  
  # Make vector to save feature set
  newProbe <- rep(NA, nFeatures[i])
  newProbe[1] <- seedProbe
  
  # Start selecting features
  for (j in 1:(nFeatures[i] - 1)){
    # Calculate correlations between seed probe and all other probes
    cor_matrix[,j] <- apply(dataMatrix_copy_fil,1,function(x){cor(x,dataMatrix_copy_fil[newProbe[j],])})
    
    # Add most uncorrelated probe to feature set
    newProbe[j+1] <- rownames(cor_matrix)[which.min(abs(cor_matrix[,j]))]
    
    # Remove all highly correlated probes: keep probes with cor < 0.9
    cor_matrix <- cor_matrix[abs(cor_matrix[,j]) < 0.5,]
    dataMatrix_copy_fil <- dataMatrix_copy_fil[rownames(cor_matrix),]
  }
  selectedProbes[[uniqueGroups[i]]] <- newProbe
}
# Record end time
t_end <- Sys.time()

# Give run time
t_end-t_start















# Make copy of data
dataMatrix_copy <- dataMatrix_M

# Make matrix to save correlations
cor_matrix <- matrix(NA, nrow(mydat_M), nFeatures-1)
rownames(cor_matrix) <- rownames(mydat_M)

# Make vector to save feature set
newProbe <- rep(NA, nFeatures)
newProbe[1] <- seedProbe

# Record start time
t_start <- Sys.time()

# Start selecting features
for (i in 1:(nFeatures - 1)){
  # Calculate correlations between seed probe and all other probes
  cor_matrix[,i] <- apply(mydat_copy,1,function(x){cor(x,mydat_copy[newProbe[i],])})
  
  # Add most uncorrelated probe to feature set
  newProbe[i+1] <- rownames(cor_matrix)[which.min(abs(cor_matrix[,i]))]
  
  # Remove all highly correlated probes: keep probes with cor < 0.9
  cor_matrix <- cor_matrix[abs(cor_matrix[,i]) < 0.5,]
  mydat_copy <- mydat_copy[rownames(cor_matrix),]
}
# Record end time
t_end <- Sys.time()

# Give run time
t_end-t_start





#*****************************************************************************#
# Probe annotation
#*****************************************************************************#

# Get annotation
#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annList <- ann@listData

probe_annotation <- data.frame(
  ID = annList$Name,
  Chr = annList$chr,
  Position = annList$pos,
  Strand = annList$strand,
  Type = annList$Type,
  CpG_MAF = annList$CpG_maf,
  Relation_to_Island = annList$Relation_to_Island,
  Gene_Name = annList$UCSC_RefGene_Name,
  Gene_ID = annList$UCSC_RefGene_Accession,
  Gene_Group = annList$UCSC_RefGene_Group,
  Regulatory_Group = annList$Regulatory_Feature_Group
)

save(probe_annotation,file = "E:/Thesis/MLData/probe_annotation.RData")





# load annotation
load("E:/Thesis/MLData/probe_annotation.RData")
library(tidyverse)
unique(probe_annotation$Regulatory_Group)



body <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "Body")]
tss200 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "TSS200")]
tss1500 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "TSS1500")]
utr5 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "5'UTR")]
utr3 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "3'UTR")]
exon1 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "1stExon")]

promotor <- unique(c(tss200,tss1500))
geneBody <- setdiff(unique(c(utr5,exon1, body, utr3)),promotor)
interGenic <- setdiff(probe_annotation$ID, unique(c(promotor, geneBody)))

probe_annotation$Class <- NA
probe_annotation$Class[probe_annotation$ID %in% promotor] <- "Promotor"
probe_annotation$Class[probe_annotation$ID %in% geneBody] <- "Gene Body"
probe_annotation$Class[probe_annotation$ID %in% interGenic] <- "Intergenic"

save(probe_annotation,file = "E:/Thesis/MLData/probe_annotation.RData")
