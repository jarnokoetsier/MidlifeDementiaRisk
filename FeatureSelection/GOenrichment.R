###############################################################################

# Perform GO enrichment analysis

###############################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("~/X_test_Cor.RData")
load("~/Data/X_test.RData")

# Load missMethyl package
library(missMethyl)

# Select CpGs
selCpGs <- rownames(X_test_Cor)
allCpGs <- rownames(X_test)

# Perform GO enrichment analysis
gst <- gometh(sig.cpg=selCpGs, all.cpg=allCpGs, collection="GO",
              array.type = "EPIC",
              plot.bias=TRUE)

# Save
save(gst, file = "GO_enrichment_Cor.RData")

# Repeat for the different feature selection methods

###############################################################################

# Plot GO enrichment results

###############################################################################

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(rrvgo)
library(ggpubr)

# Function to capatilize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation/FeatureSelection")

# Feature selection methods:
method <- c("var","varM", "varCor", "varMCor", "KS", "Cor")
methodNames <- c("Variance (\u03b2)", "Variance (M)", "Variance (\u03b2, Cor)",
                 "Variance (M, Cor)", "KS-like", "Correlation")

# Collect p-values
for (i in 1:length(method)){
  fileName <- paste0("GO_enrichment_",method[i],".RData")
  load(fileName)
  
  if (i == 1){
    temp <- gst[,c(1,2,5)]
    temp_FDR <- gst[,c(1,2,6)]
  }
  if (i > 1){
    temp <- cbind(temp, gst$P.DE)
    temp_FDR <- cbind(temp_FDR, gst$P.DE)
  }
  
}
colnames(temp) <- c("Ontology", "Term", method)
colnames(temp_FDR) <- c("Ontology", "Term", method)

# Select top GO terms
temp_fil <- NULL
for (i in 1:length(method)){
  temp_order <- temp[order(temp[,method[i]]),]
  
  temp_fil <- rbind.data.frame(temp_fil,head(temp_order,15))
}
temp_fil <- unique(temp_fil)
temp_fil_FDR <- temp_FDR[rownames(temp_fil),]

colnames(temp_fil) <- c("Ontology", "Term", methodNames)
colnames(temp_fil_FDR) <- c("Ontology", "Term", methodNames)

# 1. Biological process:

# Make similarity matrix
simMatrix_BP <- calculateSimMatrix(rownames(temp_fil)[temp_fil$Ontology == "BP"],
                                   orgdb = "org.Hs.eg.db",
                                   ont = "BP", 
                                   method = "Rel")

# reduce similar terms
reduceTerms_BP <- reduceSimMatrix(simMatrix_BP,
                                  threshold = 0.7,
                                  orgdb = "org.Hs.eg.db")

# 2. Molecular Function: 

# Make similarity matrix
simMatrix_MF <- calculateSimMatrix(rownames(temp_fil)[temp_fil$Ontology == "MF"],
                                   orgdb = "org.Hs.eg.db",
                                   ont = "MF", 
                                   method = "Rel")

# reduce similar terms
reduceTerms_MF <- reduceSimMatrix(simMatrix_MF,
                                  threshold = 0.7,
                                  orgdb = "org.Hs.eg.db")

# 3. Cellular Component:

# Make similarity matrix
simMatrix_CC <- calculateSimMatrix(rownames(temp_fil)[temp_fil$Ontology == "CC"],
                                   orgdb = "org.Hs.eg.db",
                                   ont = "CC", 
                                   method = "Rel")

# reduce similar terms
reduceTerms_CC <- reduceSimMatrix(simMatrix_CC,
                                  threshold = 0.7,
                                  orgdb = "org.Hs.eg.db")


# Collect reduced terms
reduceTerms <- rbind.data.frame(reduceTerms_BP, reduceTerms_MF, reduceTerms_CC)

temp_fil <- left_join(temp_fil, reduceTerms[,c("go","term", "parent")], by = c("Term" = "term"))
temp_fil_FDR <- left_join(temp_fil_FDR, reduceTerms[,c("go","term", "parent")], by = c("Term" = "term"))
temp_fil$parent[is.na(temp_fil$parent)] <- temp_fil$parent[2]

FDR <- gather(temp_fil_FDR[,3:8])
rownames(temp_fil) <- paste0(temp_fil$Term, " (", temp_fil$Ontology, ")")
rownames(temp_fil) <- firstup(rownames(temp_fil))

# prepare data for plotting
plotDF <- gather(temp_fil[,3:8])
plotDF$Ontology <- rep(temp_fil$Ontology, 6)
plotDF$Term <- rep(temp_fil$Term, 6)
plotDF$FDR <- FDR$value
plotDF$Sig <- ifelse(plotDF$FDR < 0.05, "Yes", "No")
plotDF$Name <- paste0(plotDF$Term, " (", plotDF$Ontology, ")")
plotDF$Name <- firstup(plotDF$Name)
plotDF$Parent <- rep(temp_fil$parent,6)

plotDF$key <- factor(plotDF$key, levels = methodNames)

# Perform hierarchical clustering on the GO terms
clusters <- hclust(dist(-log(temp_fil[,3:8])), method = "ward.D2")
order <- clusters$labels[clusters$order]
plotDF$Name <- factor(plotDF$Name, levels = order)

# Main plot
main <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Name, fill = -log(value), color = Sig),
            width = 0.9, height = 0.9, linewidth = 0.5) +
  scale_fill_viridis_c() +
  scale_color_manual(values = c("white","black")) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = "-log p-value") +
  guides(color = "none") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

# Side colors: GO term similarity
colors <- rev(RColorBrewer::brewer.pal(n = 11, "Set3"))
side <- ggplot(unique(plotDF)) +
  geom_tile(aes(x = 1, y = Name, fill = Parent),
            width = 0.9, height = 0.9, linewidth = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = c("white","black")) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# Combine plots
p <- ggarrange(side, main, ncol = 2, nrow = 1, align = "h", widths = c(5.5,5))

# Save figure
ggsave(p, file = "heatmap_GOenrichment_featureselection.png", width = 9, height = 10)


