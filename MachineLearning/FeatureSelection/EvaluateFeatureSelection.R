# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Load S-scored data
files <- list.files('X_S')
for (f in files){
  load(paste0("X_S/",f))
}
# Load var data
files <- list.files('X_var')
for (f in files){
  load(paste0("X_var/",f))
}

# Load packages
library(tidyverse)


###############################################################################

# 1. S-score

###############################################################################



#*****************************************************************************#
# Beta vs sd
#*****************************************************************************#
Mvalues_S <- log2(X_nonTest_S/(1 - X_nonTest_S))
Mvalues_var <- log2(X_nonTest_var/(1 - X_nonTest_var))

plot_S <- data.frame(
  meanBeta = rowMeans(X_nonTest_S),
  sdBeta = apply(X_nonTest_S, 1, sd),
  meanM  = rowMeans(Mvalues_S),
  sdM = apply(Mvalues_S, 1, sd)
)

plotVar <- data.frame(
  meanBeta = rowMeans(X_nonTest_var),
  sdBeta = apply(X_nonTest_var, 1, sd),
  meanM  = rowMeans(Mvalues_var),
  sdM = apply(Mvalues_var, 1, sd)
)

ggplot() +
  geom_point(data = plot_S, aes(x = meanBeta, y = sdBeta, color = "S-score")) +
  geom_point(data = plotVar, aes(x = meanBeta, y = sdBeta, color = "Variance")) +
  scale_color_brewer(palette = "Dark2")

ggplot() +
  geom_point(data = plot_S, aes(x = meanM, y = sdM, color = "S-score")) +
  geom_point(data = plotVar, aes(x = meanM, y = sdM, color = "Variance")) +
  scale_color_brewer(palette = "Dark2")
#*****************************************************************************#
# PCA
#*****************************************************************************#

# Perform PCA
pcaList <-  prcomp(t(X_nonTest_S),        
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
  ggtitle("S-score Selection") +
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

# Load S-scored data
files <- list.files('X_var')
for (f in files){
  load(paste0("X_var/",f))
}
#*****************************************************************************#
# PCA
#*****************************************************************************#

# Perform PCA
pcaList_var <-  prcomp(t(X_nonTest_var),        
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
PCAscores_var$Sex <- ifelse(PCAscores_var$Sex == 1, "Male", "Female")
p_12 <- ggplot(data = PCAscores_var, aes(x = PC1, y = PC2)) +
  stat_ellipse(geom = "polygon",
               aes(fill = factor(Sex)),
               type = "norm", 
               alpha = 0.25,
               level = 0.99) +
  geom_point(alpha = 0.9, size = 2, aes(color = Neu)) +
  xlab(paste0("PC1 (", explVar_var[1],"%)")) +
  ylab(paste0("PC2 (", explVar_var[2],"%)")) +
  ggtitle("Variance Selection") +
  labs(color = "Neutrophils", fill = "Sex") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text()) +
  scale_color_viridis_c() +
  scale_fill_brewer(palette = "Set1")

# save plot
ggsave(p_12, file = "PCAplot_1vs2_var.png", width = 8, height = 6)

#=============================================================================#
# PCA correlations
#=============================================================================#

PCs <- PCAscores_var[,1:5]
colnames(PCs) <- paste0(colnames(PCs), " (", explVar_var[1:5], "%)")
meta <- PCAscores_var[,c("Age", "Sex", colnames(cellType))]
meta <- meta[,-c(9,10)]
meta$Sex <- ifelse(meta$Sex == "Male", 0,1)
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
  ggtitle("Variance Selection") +
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


