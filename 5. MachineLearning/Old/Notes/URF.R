
install.packages("randomForest")


library(randomForest)
library(tidyverse)
library(RColorBrewer)
library(heatmaply)
library(corrr)
setwd("E:/Thesis/EXTEND/Methylation")

#******************************************************************************#
# 2.1. Unsupervised random forest
#******************************************************************************#
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

cellType_scaled <- t((t(cellType[,1:6]) - colMeans(cellType[,1:6]))/(apply(cellType[,1:6],2,sd)))
diss <- as.matrix(dist(cellType_scaled, method = "euclidean"))



# Center data (double mean centering)
diss_centered <- diss - rowMeans(diss) - colMeans(diss) + mean(diss)

# Perform PCA
pcoaList <-  prcomp(diss_centered,        
                    retx = TRUE,
                    center = FALSE,
                    scale = FALSE)


# Explained variance
explVarPCoA <- round(((pcoaList$sdev^2)/sum(pcoaList$sdev^2))*100,2)

# Get PCoA scores
plotPCoA_scores <- as.data.frame(pcoaList$x[,1:10])
plotPCoA_scores$ID <- rownames(plotPCoA_scores)

# Combine PCoA scores with cell type information
plotPCoA_scores <- inner_join(plotPCoA_scores, cellType, by = c("ID" = "ID"))

# Combine PCoA scores with cluster information
plotPCoA_scores <- inner_join(plotPCoA_scores, dat, by = c("ID" = "Basename"))

# Make PCoA score plots:

# 1) Color by URF-based clusters
PCoA <- ggplot(plotPCoA_scores, aes(x = PC1, y = PC3, color = CD4T)) +
  geom_point(size = 2, alpha = 0.9) +
  xlab(paste0("PCo1 (", explVarPCoA[1],"%)")) +
  ylab(paste0("PCo2 (", explVarPCoA[2],"%)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_color_viridis_c()


PCs <- plotPCoA_scores[,1:7]
meta <- plotPCoA_scores[,12:17]
meta <- meta[,sapply(meta, is.numeric)]


corDF <- expand.grid(colnames(PCs), colnames(meta))
colnames(corDF) <- c("PCs", "Meta")
corDF$Correlation <- rep(NA, nrow(corDF))
for (i in 1:nrow(corDF)){
  corDF$Correlation[i] <- cor(PCs[,corDF$PCs[i]], meta[,corDF$Meta[i]], method = "spearman")
}

p <- ggplot() +
  geom_point(data = corDF, aes(x = PCs, y = Meta, color = Correlation, size = abs(Correlation))) +
  labs(color = "Spearman\nCorrelation", size = "|Spearman\nCorrelation|") +
  scale_color_viridis_c() +
  theme(axis.title = element_blank())

ggsave(p,file = "correlationPlot.png", width = 8, height = 6)


