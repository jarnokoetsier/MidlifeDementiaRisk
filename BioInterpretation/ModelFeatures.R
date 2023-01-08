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
library(clusterProfiler)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load cell type composition and phenotype data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Load model information
Score = "CAIDE1"
FeatureSelection = "Cor"
load(paste0("CV_CAIDE1/CV_", Score, "_", FeatureSelection,".RData"))



coefs_finalModel <- as.data.frame(as.matrix(coef(finalModel)))
probes <- rownames(coefs_finalModel)[coefs_finalModel$s0!=0][-1]

test <- coefs_finalModel$s0
names(test) <- rownames(coefs_finalModel)
probes <- names(tail(sort(test),101)[-101])

# Get gene annotations for the probes
load("Probe2Gene_final_ann.RData")
Probe2Gene_all <- unique(Probe2Gene_final[,c(2,4)])



#read GMT file
gmt <- clusterProfiler::read.gmt.wp("wikipathways-20220810-gmt-Homo_sapiens.gmt.txt")

# link pathway to gene
path2gene <- gmt[,c("wpid", "gene")]
path2name <- gmt[,c("wpid", "name")]

geneSet <- Probe2Gene_final$EntrezID[Probe2Gene_final$CpG %in% probes]

WP <- enricher(gene = as.character(unique(geneSet)),
               universe = as.character(unique(Probe2Gene_final$EntrezID)),
               pAdjustMethod = "fdr",
               pvalueCutoff = 1,
               TERM2GENE = path2gene,
               TERM2NAME = path2name)

result_WP_CAIDE1 <- WP@result
result_WP_CAIDE1_copy <- WP@result[,c(1,2,5)]

# Number of permutations
nPerm <- 100

# All CpGs
all_cpg <- unique(Probe2Gene_final$CpG)

# Number of cpgs to select during each permutation
n_cpg <- length(probes)

#set.seed(59807)
for (i in 1:nPerm){
  
  # Get random miRNAs
  cpg <- sample(all_cpg, n_cpg)
  
  # Get asssociated genes
  geneSet_perm <- as.character(Probe2Gene_final$EntrezID[Probe2Gene_final$CpG %in% cpg])
  
  # Perform WP ORA
  WP <- enricher(gene = geneSet_perm,
                 universe = as.character(unique(Probe2Gene_final$EntrezID)),
                 pAdjustMethod = "fdr",
                 pvalueCutoff = Inf,
                 qvalueCutoff = Inf,
                 TERM2GENE = path2gene,
                 TERM2NAME = path2name)
  
  # Get results
  results <- WP@result[,c(1,5)]
  colnames(results) <- c("ID", paste0("Perm",i))
  
  # Combine results
  result_WP_CAIDE1_copy <- left_join(result_WP_CAIDE1_copy, results, by = c("ID" = "ID"))
}


save(result_CAIDE1_early_copy, file = "results_WP_CAIDE1_copy.RData")





# Get permutation results
perm_results <- -log10(result_WP_CAIDE1_copy[,4:(nPerm+3)])
perm_results[is.na(perm_results)] <- 0

# Count number of times with more enrichment
test <- rowSums(perm_results > -log10(result_WP_CAIDE1_copy[,3]))/nPerm

# Combine into data.frame
Output <- cbind.data.frame(result_WP_CAIDE1_copy[,1:3], test)

# Calculate FDR
Output$FDR <- p.adjust(test, method = "fdr")

# Change column names
colnames(Output) <- c("ID", "Description", "pvalue (no perm)", "pvalue (perm)", "FDR")

# Order the rows by pvalue
Output <- arrange(Output, by = `pvalue (perm)`)

# Write output
save(perm_results, file = "perm_results_Early.RData")
write.csv(Output, file = "WP_Early.csv")

