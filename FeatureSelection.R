
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

# Convert to M-values
mydat1_fil <- mydat1[rowSums((mydat1 > 0) & (mydat1 < 1)) == 919, ]
mydat_M <- log2(mydat1_fil/(1 + mydat1_fil))

# Number of features to select
nFeatures <- 500

# Seed probe
seedProbe <- rownames(mydat_M)[323]

# Make copy of data
mydat_copy <- mydat_M

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

test <- data.frame(matrix(unlist(ann@listData), nrow=length(ann@listData), byrow=TRUE),stringsAsFactors=FALSE)

test <- do.call(rbind.data.frame, ann@listData)
test <- lapply(ann@listData, cbind.data.frame)

annDF <- unlist(ann@listData)

flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
