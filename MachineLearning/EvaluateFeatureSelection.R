# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load data
load("dataMatrix_S.RData")
load("dataMatrix_var.RData")
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Load packages
library(tidyverse)


###############################################################################

# 1. S-score

###############################################################################

#*****************************************************************************#
# PCA
#*****************************************************************************#

# Perform PCA
pcaList <-  prcomp(t(dataMatrix_S),        
                   retx = TRUE,
                   center =TRUE,
                   scale = TRUE,
                   rank. = 10)

# Save pcaList object
save(pcaList, file = "pcaList_S.RData")
load("pcaList_S.RData")

# Get PCA scores
PCAscores <- as.data.frame(pcaList$x)
PCAscores$ID <- rownames(PCAscores)
PCAscores <- inner_join(PCAscores, dat, by = c("ID" = "Basename"))

# Combine with cell type composition
PCAscores <- inner_join(PCAscores,cellType, by = c("ID" = "ID"))

# Get explained variance
explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)

#=============================================================================#
# PCA score plot
#=============================================================================#

# Make PCA score: PC1 vs PC2
p_12 <- ggplot(data = PCAscores, aes(x = PC1, y = PC2)) +
  stat_ellipse(geom = "polygon",
               fill = "red",
               type = "norm", 
               alpha = 0.25,
               level = 0.95) +
  stat_ellipse(geom = "polygon",
               color = "red",
               alpha = 0,
               linetype = 2,
               type = "norm",
               level = 0.99) +
  geom_point(alpha = 0.9, size = 2, aes(color = Sex)) +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  labs(color = "CD8 T-cells") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text()) +
  scale_color_viridis_c()

# save plot
ggsave(p_12, file = "PCAplot_1vs2.png", width = 8, height = 6)

#=============================================================================#
# PCA correlations
#=============================================================================#

PCs <- PCAscores[,1:5]
colnames(PCs) <- paste0(colnames(PCs), " (", explVar[1:5], "%)")
meta <- PCAscores[,c("Age", "Sex", colnames(cellType))]
meta <- meta[,-c(9,10)]
#meta$Sex <- ifelse(meta$Sex == "Male", 0,1)
colnames(meta) <- c("Age", "Sex", "CD8 T-cells", "CD4 T-cells", "NK cells", "B-cells", "Monocytes",
                    "Neutrophils")

corDF <- expand.grid(colnames(PCs), colnames(meta))
colnames(corDF) <- c("PCs", "Meta")
corDF$Correlation <- rep(NA, nrow(corDF))
for (i in 1:nrow(corDF)){
  corDF$Correlation[i] <- cor(PCs[,corDF$PCs[i]], meta[,corDF$Meta[i]], method = "spearman")
}

p <- ggplot() +
  geom_point(data = corDF, aes(x = PCs, y = Meta, color = Correlation, size = abs(Correlation))) +
  labs(color = "Spearman\nCorrelation", size = "|Spearman\nCorrelation|") +
  guides(size = "none") +
  ggtitle("S-score") +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                   face = "bold",
                                   size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                      size = 10))

ggsave(p,file = "correlationPlot_S.png", width = 8, height = 6)


###############################################################################

# 2. Variance

###############################################################################

#*****************************************************************************#
# PCA
#*****************************************************************************#

# Perform PCA
pcaList_var <-  prcomp(t(dataMatrix_var),        
                   retx = TRUE,
                   center =TRUE,
                   scale = TRUE,
                   rank. = 10)

# Save pcaList object
save(pcaList_var, file = "pcaList_var.RData")

# Get PCA scores
PCAscores_var <- as.data.frame(pcaList_var$x)
PCAscores_var$ID <- rownames(PCAscores_var)
PCAscores_var <- inner_join(PCAscores_var, dat, by = c("ID" = "Basename"))

# Combine with cell type composition
PCAscores_var <- inner_join(PCAscores_var,cellType, by = c("ID" = "ID"))

# Get explained variance
explVar_var <- round(((pcaList_var$sdev^2)/sum(pcaList_var$sdev^2))*100,1)

#=============================================================================#
# PCA score plot
#=============================================================================#

# Make PCA score: PC1 vs PC2
p_12 <- ggplot(data = PCAscores_var, aes(x = PC1, y = PC2)) +
  stat_ellipse(geom = "polygon",
               fill = "red",
               type = "norm", 
               alpha = 0.25,
               level = 0.95) +
  stat_ellipse(geom = "polygon",
               color = "red",
               alpha = 0,
               linetype = 2,
               type = "norm",
               level = 0.99) +
  geom_point(alpha = 0.9, size = 2, aes(color = Sex)) +
  xlab(paste0("PC1 (", explVar_var[1],"%)")) +
  ylab(paste0("PC2 (", explVar_var[2],"%)")) +
  labs(color = "CD8 T-cells") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text()) +
  scale_color_viridis_c()

# save plot
ggsave(p_12, file = "PCAplot_1vs2.png", width = 8, height = 6)

#=============================================================================#
# PCA correlations
#=============================================================================#

PCs <- PCAscores_var[,1:5]
colnames(PCs) <- paste0(colnames(PCs), " (", explVar_var[1:5], "%)")
meta <- PCAscores_var[,c("Age", "Sex", colnames(cellType))]
meta <- meta[,-c(9,10)]
#meta$Sex <- ifelse(meta$Sex == "Male", 0,1)
colnames(meta) <- c("Age", "Sex", "CD8 T-cells", "CD4 T-cells", "NK cells", "B-cells", "Monocytes",
                    "Neutrophils")

corDF <- expand.grid(colnames(PCs), colnames(meta))
colnames(corDF) <- c("PCs", "Meta")
corDF$Correlation <- rep(NA, nrow(corDF))
for (i in 1:nrow(corDF)){
  corDF$Correlation[i] <- cor(PCs[,corDF$PCs[i]], meta[,corDF$Meta[i]], method = "spearman")
}

p <- ggplot() +
  geom_point(data = corDF, aes(x = PCs, y = Meta, color = Correlation, size = abs(Correlation))) +
  labs(color = "Spearman\nCorrelation", size = "|Spearman\nCorrelation|") +
  guides(size = "none") +
  ggtitle("Variance") +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))

ggsave(p,file = "correlationPlot_Var.png", width = 8, height = 6)


