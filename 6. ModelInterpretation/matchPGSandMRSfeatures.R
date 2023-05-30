
# Load packages
library(caret)
library(GenomicRanges)
library(tidyverse)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("~/allModels.RData")
load("~/Data/Probe2Gene_final_ann.RData")

#*****************************************************************************#
#   Common genes
#*****************************************************************************#

# Collect gene location information
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
broads <- GenomicFeatures::genes(txdb)

# Get top 500 features of MRS model
features_MRS <- varImp(allModels$Education, scale = FALSE)$importance
features_MRS <- rownames(head(arrange(features_MRS, -1*Overall),500))

# Get associated genes for each CpG
genes_MRS <- unique(Probe2Gene_final$EntrezID[Probe2Gene_final$CpG %in% features_MRS])

# Get top 500 features of PGS model
features_PGS <- read.delim("Data/bayesr.effects", sep = " ")
features_PGS$Chr <-  unlist(lapply(str_split(features_PGS$Predictor, ":"), `[[`, 1))
features_PGS$Position <- unlist(lapply(str_split(features_PGS$Predictor, ":"), `[[`, 2))
features_PGS1 <- head(arrange(features_PGS, -1*abs(Model58)),500)

# Get nearest gene for each SNP
x <- GRanges(
  seqnames = paste0("chr",features_PGS1$Chr),
  ranges = IRanges(features_PGS1$Position, names = features_PGS1$Predictor),
  strand = rep("*", length(features_PGS1$Chr)))

nearestGenes <- nearest(x,broads)
nearestGenes <- broads[nearestGenes,]@ranges@NAMES

# Overlap between genes associated with CpGs and SNPs
overlap <- intersect(genes_MRS, nearestGenes)

# Annotate ENTREZ ID with gene symbol
library(org.Hs.eg.db)
geneSymbol <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = overlap,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "ENTREZID")
# Get unique symbols only
geneSymbol <- unique(geneSymbol)
save(geneSymbol, file = "OverlapGenes_AE.RData")


# Perform permutation to test for significance
nGenes_all <- NULL
set.seed(123)
for (i in 1:10000){
  
  # Select 500 random SNPs and their associated genes
  PGS_rand <- features_PGS[sample(1:nrow(features_PGS1),500),]
  x <- GRanges(
    seqnames = paste0("chr",PGS_rand$Chr),
    ranges = IRanges(PGS_rand$Position, names = PGS_rand$Predictor),
    strand = rep("*", length(PGS_rand$Chr)))
  
  nearestGenes <- nearest(x,broads)
  nearestGenes <- broads[nearestGenes,]@ranges@NAMES
  
  # Select 500 random CpGs and their associated genes
  all_CpG <- unique(Probe2Gene_final$CpG)
  MRS_rand <- all_CpG[sample(1:length(all_CpG),500)]
  genes_MRS <- unique(Probe2Gene_final$EntrezID[Probe2Gene_final$CpG %in% MRS_rand])
  
  # Determine overlap
  temp <- length(intersect(genes_MRS, nearestGenes))
  nGenes_all <- c(nGenes_all,temp)
}
# Save permutation results
save(nGenes_all, file = "nGenes_all.RData")

# Determine significance (more than or equal to 7 overlapping genes)
sum(nGenes_all >= 7)

# Plot histogram of permutation results
plotDF <- data.frame(N = nGenes_all)
p <- ggplot(plotDF) +
  geom_density(aes(x = N, y = after_stat(density)), adjust = 5, 
               fill = "grey", alpha = 0.3, color = "white") +
  geom_histogram(aes(x = N, y = after_stat(density)*0.5), fill = "grey", color = "black", bins = 15) +
  xlab("Number of common genes") +
  ylab("Count") +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3), labels = c(0,1000,2000,3000)) +
  theme_classic()

# Save plot
ggsave(p, file = "overlapGenes_MRSPGS.png", width = 8, height = 5)


#*****************************************************************************#
#   mQTLs
#*****************************************************************************#

# Load meQTLs (ARIES mQTL database)
test <- read.csv("meQTLs.csv")

# For all PGSs, check if they are in the top 500 SNPs and whether they are meQTLs of
# the top 500 CpGs of MRS model
stat <- data.frame(Predictor = features_PGS$Predictor,
                   Top = rep("No", length(features_PGS$Predictor)),
                   meQTL = rep("No", length(features_PGS$Predictor)))

stat$meQTL[stat$Predictor %in% paste0(test$SNP.Chr, ":", test$SNP.Pos)] <- "Yes"
stat$Top[stat$Predictor %in% features_PGS1$Predictor] <- "Yes"

# there is a low expected value (< 5): perform Fisher's exact test
chisq.test(stat$meQTL, stat$Top)$expected
fisher.test(stat$meQTL, stat$Top)
table(stat$meQTL, stat$Top)