library(tidyverse)
gmt <- clusterProfiler::read.gmt.wp("E:/RTTproject/ExosomeAnalysis/4. Pathway Analysis/wikipathways-20220810-gmt-Homo_sapiens.gmt.txt")
load("E:/Thesis/EXTEND/OverlapGenes_AE.RData")

geneSymbol <- c("TCF4", "SORCS1", "SEMA6D", "PYM1", "PNMA8A", "PHF2", "PHACTR3",
                "MECOM", "LRP10", "KRT82", "EMP1", "DPYD", "DHFR", "CORO6", "CDH4",
                "BAG5")


test <- read.csv("C:/Users/Gebruiker/Desktop/H. sapiens (3) default node.csv")
geneSymbol <- test$gene.name
library(org.Hs.eg.db)
geneSymbol <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = geneSymbol,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
# Get unique symbols only
library(clusterProfiler)
ego <- enrichGO(gene          = geneSymbol$ENTREZID,
                universe      = NULL,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

test <- ego@result


wp_fil <- inner_join(gmt, geneSymbol, by = c("gene" = "ENTREZID"))

nodes <- data.frame(id = c(unique(wp_fil$name), unique(wp_fil$SYMBOL)),
                    Group = c(rep("Pathway", length(unique(wp_fil$name))),
                              rep("Gene", length(unique(wp_fil$SYMBOL))))
)

edges <- data.frame(source = wp_fil$name,
                    target  = wp_fil$SYMBOL)

RCy3::createNetworkFromDataFrames(nodes= nodes, edges = edges, title="Test", collection="Test")

