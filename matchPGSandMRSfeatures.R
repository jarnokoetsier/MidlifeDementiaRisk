library(caret)
library(GenomicRanges)
library(tidyverse)

# Clear workspace and console
rm(list = ls())
cat("\014") 

load("~/allModels.RData")
load("~/Data/Probe2Gene_final_ann.RData")


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
broads <- GenomicFeatures::genes(txdb)


features_MRS <- varImp(allModels$Education, scale = FALSE)$importance
features_MRS <- rownames(head(arrange(features_MRS, -1*Overall),500))
genes_MRS <- unique(Probe2Gene_final$EntrezID[Probe2Gene_final$CpG %in% features_MRS])

features_PGS <- read.delim("Data/bayesr.effects", sep = " ")
features_PGS$Chr <-  unlist(lapply(str_split(features_PGS$Predictor, ":"), `[[`, 1))
features_PGS$Position <- unlist(lapply(str_split(features_PGS$Predictor, ":"), `[[`, 2))
features_PGS1 <- head(arrange(features_PGS, -1*abs(Model58)),500)


x <- GRanges(
  seqnames = paste0("chr",features_PGS1$Chr),
  ranges = IRanges(features_PGS1$Position, names = features_PGS1$Predictor),
  strand = rep("*", length(features_PGS1$Chr)))

nearestGenes <- nearest(x,broads)
nearestGenes <- broads[nearestGenes,]@ranges@NAMES


overlap <- intersect(genes_MRS, nearestGenes)

features_PGS1$Gene <- nearestGenes
PGS_DF <- features_PGS1[which(nearestGenes %in% overlap),]


MRS_DF <- Probe2Gene_final[(Probe2Gene_final$EntrezID %in% overlap) & (Probe2Gene_final$CpG %in% features_MRS),]
MRS_DF <- unique(MRS_DF[,c(1,2,4)])
coefs <- as.matrix(coef(allModels$Education$finalModel,allModels$Education$finalModel$lambdaOpt))
MRS_DF$Coef <- coefs[MRS_DF$CpG,1]
sum(coefs[,1] != 0)


total_DF <- inner_join(PGS_DF, MRS_DF, by = c("Gene" = "EntrezID"))

# Annotate ENTREZ ID with gene symbol
library(org.Hs.eg.db)
geneSymbol <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = overlap,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "ENTREZID")
# Get unique symbols only
geneSymbol <- unique(geneSymbol)

save(geneSymbol, file = "OverlapGenes_AE.RData")


# Perform permutation

nGenes_all <- NULL
set.seed(123)
for (i in 1:10000){
  PGS_rand <- features_PGS[sample(1:nrow(features_PGS1),500),]
  x <- GRanges(
    seqnames = paste0("chr",PGS_rand$Chr),
    ranges = IRanges(PGS_rand$Position, names = PGS_rand$Predictor),
    strand = rep("*", length(PGS_rand$Chr)))
  
  nearestGenes <- nearest(x,broads)
  nearestGenes <- broads[nearestGenes,]@ranges@NAMES
  
  all_CpG <- unique(Probe2Gene_final$CpG)
  MRS_rand <- all_CpG[sample(1:length(all_CpG),500)]
  genes_MRS <- unique(Probe2Gene_final$EntrezID[Probe2Gene_final$CpG %in% MRS_rand])
  
  
  temp <- length(intersect(genes_MRS, nearestGenes))
  
  nGenes_all <- c(nGenes_all,temp)
}
save(nGenes_all, file = "nGenes_all.RData")

sum(nGenes_all >= 7)


plotDF <- data.frame(N = nGenes_all)

p <- ggplot(plotDF) +
  geom_density(aes(x = N, y = after_stat(density)), adjust = 5, 
               fill = "grey", alpha = 0.3, color = "white") +
  geom_histogram(aes(x = N, y = after_stat(density)*0.5), fill = "grey", color = "black", bins = 15) +
  xlab("Number of common genes") +
  ylab("Count") +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3), labels = c(0,1000,2000,3000)) +
  theme_classic()

ggsave(p, file = "overlapGenes_MRSPGS.png", width = 8, height = 5)
