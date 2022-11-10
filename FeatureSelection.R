
# Set working directory
setwd("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA")

# Load normalized training data
load("E:/Thesis/mydat1.RData")
load("E:/Thesis/pheno.RData")

# Convert to M-values
mydat1_fil <- mydat1[rowSums((mydat1 > 0) & (mydat1 < 1)) == 919, ]
mydat_M <- log2(mydat1_fil/(1 + mydat1_fil))

# Number of features to select
nFeatures <- 10

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
  cor_matrix <- cor_matrix[abs(cor_matrix[,i]) < 0.9,]
  mydat_copy <- mydat_copy[rownames(cor_matrix),]
}
# Record end time
t_end <- Sys.time()

# Give run time
t_end-t_start


# Get annotation
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

annDF <- unlist(ann@listData)

flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
