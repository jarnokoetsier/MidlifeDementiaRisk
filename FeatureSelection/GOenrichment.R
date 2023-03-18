# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation/FeatureSelection")

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

FDR <- gather(temp_fil_FDR[,3:ncol(temp_fil_FDR)])
rownames(temp_fil) <- paste0(temp_fil$Term, " (", temp_fil$Ontology, ")")
rownames(temp_fil) <- firstup(rownames(temp_fil))

plotDF <- gather(temp_fil[,3:ncol(temp_fil)])
plotDF$Ontology <- rep(temp_fil$Ontology, ncol(temp_fil) - 2)
plotDF$Term <- rep(temp_fil$Term, ncol(temp_fil) - 2)
plotDF$FDR <- FDR$value
plotDF$Sig <- ifelse(plotDF$FDR < 0.05, "Yes", "No")
plotDF$Name <- paste0(plotDF$Term, " (", plotDF$Ontology, ")")
plotDF$Name <- firstup(plotDF$Name)

plotDF$key <- factor(plotDF$key, levels = methodNames)


clusters <- hclust(dist(-log(temp_fil[,3:ncol(temp_fil)])), method = "ward.D2")
order <- clusters$labels[clusters$order]
plotDF$Name <- factor(plotDF$Name, levels = order)

p <- ggplot(plotDF) +
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
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

ggsave(p, file = "heatmap_GOenrichment_featureselection.png", width = 8, height = 10)



fileName <- paste0("GO_enrichment_","S",".RData")
load(fileName)
