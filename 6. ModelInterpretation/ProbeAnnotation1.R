# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(foreach)
library(doParallel)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load probe annotation data
load("probe_annotation.RData")
load("Probe2Gene_final_ann.RData")

###############################################################################

# Make probe annotation

###############################################################################

# Filter for probes
probe_annotation <- probe_annotation[(probe_annotation$Chr != "chrY") & (probe_annotation$Chr != "ChrX"),]


# Get associated gene(s)
n <- nrow(probe_annotation)

# Register cores for parallel computing
# start from 500001
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)
Probe2Gene <- foreach (i = 700001:n, 
                       .packages = "tidyverse", 
                       .combine = rbind.data.frame, 
                       .inorder = FALSE) %dopar% {
  temp <- data.frame(
    Gene = unlist(str_split(probe_annotation$Gene_Name[i], ";")),
    CpG = rep(probe_annotation$ID[i],length(unlist(str_split(probe_annotation$Gene_Name[i], ";")))),
    Association = unlist(str_split(probe_annotation$Gene_Group[i], ";"))
  )
  return(temp)
}
# Stop clusters
stopCluster(cl)

Probe2Gene_final <- rbind.data.frame(Probe2Gene_final,Probe2Gene)
save(Probe2Gene_final, file = "Probe2Gene_final.RData")


###############################################################################

# Analyse probe annotation

###############################################################################

library(GenomicRanges)

#=============================================================================#
# Get intergenic CpGs
intergenic <- Probe2Gene_final$CpG[Probe2Gene_final$Gene == ""]
intergenic <- probe_annotation[probe_annotation$ID %in% intergenic,c(1,2,3,4)]

# Get nearest genes
x <- GRanges(
  seqnames = intergenic$Chr,
  ranges = IRanges(intergenic$Position, names = intergenic$ID),
  strand = intergenic$Strand)

nearestGenes <- nearest(x)
nearestGenes <- broads[nearestGenes,]@ranges@NAMES

# Annotate ENTREZ ID with gene symbol
library(org.Hs.eg.db)
geneSymbol <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = nearestGenes,
                                  columns = c("ENTREZID", "SYMBOL"),
                                  keytype = "ENTREZID")
# Get unique symbols only
geneSymbol <- unique(geneSymbol)

# get gene symbol for each gene
rownames(geneSymbol) <- geneSymbol$ENTREZID
geneSymbol <- geneSymbol[nearestGenes,]
geneSymbol$SYMBOL[is.na(geneSymbol$SYMBOL)] <- geneSymbol$ENTREZID[is.na(geneSymbol$SYMBOL)]

# Check whether genes are in correct order
all(nearestGenes == geneSymbol$ENTREZID)
nearestGenes <- geneSymbol$SYMBOL 

# Combine into dataframe
Probe2Gene_intergenic <- data.frame(Gene = nearestGenes,
                                    CpG = intergenic$ID,
                                    Association = rep("Intergenic", length(intergenic$ID)))

Probe2Gene_nonIntergenic <- Probe2Gene_final[Probe2Gene_final$Gene != "",]


Probe2Gene_final <- rbind.data.frame(Probe2Gene_nonIntergenic, Probe2Gene_intergenic)
save(Probe2Gene_final, file = "Probe2Gene_final_ann.RData")
#==============================================================================#

# Get ENTREZ IDs
library(org.Hs.eg.db)
EntrezID <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = unique(Probe2Gene_final$Gene),
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")



id <- rep(NA, length(Probe2Gene_final$Gene))
for (gene in 1:length(Probe2Gene_final$Gene)){
  temp <- EntrezID$ENTREZID[EntrezID$SYMBOL == Probe2Gene_final$Gene[gene]]
  id[gene] <- temp[1]
}

Probe2Gene_final$EntrezID <- id
Probe2Gene_final$EntrezID[is.na(Probe2Gene_final$EntrezID)] <- Probe2Gene_final$Gene[is.na(Probe2Gene_final$EntrezID)]



###############################################################################

# Make plots

###############################################################################


Probe2Gene_all <- unique(Probe2Gene_final[,1:2])

Probe2Gene_Promotor <- unique(Probe2Gene_final[(Probe2Gene_final$Association == "TSS200") |
                                                 (Probe2Gene_final$Association == "TSS1500"),1:2])

Probe2Gene_GeneBody <- unique(Probe2Gene_final[(Probe2Gene_final$Association == "Body") |
                                                 (Probe2Gene_final$Association == "1stExon") |
                                                 (Probe2Gene_final$Association == "ExonBnd") |
                                                 (Probe2Gene_final$Association == "5'UTR") |
                                                 (Probe2Gene_final$Association == "3'UTR"),1:2])

Probe2Gene_Intergenic <- unique(Probe2Gene_final[(Probe2Gene_final$Association == "Intergenic"),1:2])

tail(sort(table(Probe2Gene_all$Gene)))


p <- ggplot()+
  geom_histogram(data = as.data.frame(table(Probe2Gene_all$Gene)),
                 aes(x = Freq, fill = "Total"), 
                 bins = 50, alpha = 0.5, color = "white") +
  geom_histogram(data = as.data.frame(table(Probe2Gene_Promotor$Gene)),
                 aes(x = Freq, fill = "Promotor"), 
                 bins = 50, alpha = 0.5, color = "white") +
  geom_histogram(data = as.data.frame(table(Probe2Gene_GeneBody$Gene)),
                 aes(x = Freq,  fill = "Gene Body"), 
                 bins = 50, alpha = 0.5, color = "white") +
  geom_histogram(data = as.data.frame(table(Probe2Gene_Intergenic$Gene)),
                 aes(x = Freq,  fill = "Intergenic"), 
                 bins = 50, alpha = 0.5, color = "white") +
  xlim(0,100)+
  xlab("# CpGs per gene") +
  ylab("Count") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave(p, file = "nCpGsPerGene.png", width = 8, height = 6)
  
ggplot()+
  geom_histogram(data = as.data.frame(table(Probe2Gene_all$Gene)),
                 aes(x = Freq), 
                 bins = 50, alpha = 0.5, color = "white") +
  xlim(0,100)+
  xlab("# CpGs per gene") +
  ylab("Count") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.title = element_blank())





