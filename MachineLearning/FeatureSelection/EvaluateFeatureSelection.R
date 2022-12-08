# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Cell type composition data
load("cellType.RData")

# Phenotype data
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Meta data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Feature selection data
FeatureSelection = "S"
files <- list.files(paste0("X", FeatureSelection))
for (f in files){
  load(paste0("X", FeatureSelection, "/",f))
}


# Make subtitle of figure
if (FeatureSelection == "Non"){
  subtitle = "No Feature Selection"
}
if (FeatureSelection == "var"){
  subtitle = "Variance (\u03b2)-based Feature Selection"
}
if (FeatureSelection == "varM"){
  subtitle = "Variance (M)-based Feature Selection"
}
if (FeatureSelection == "S"){
  subtitle = "S-score-based Feature Selection"
}

###############################################################################

# 1. Mean Beta-value vs SD Beta-value

###############################################################################

# Load all beta/M and their sd
load("plot_all.RData")
plot_all$ID <- rownames(plot_all)

# Get beta- and M-values
Bvalues <- X_nonTest_S
Mvalues <- log2(Bvalues/(1 - Bvalues))

# Format data
plotDF <- data.frame(
  ID = rownames(Bvalues),
  meanBeta = rowMeans(Bvalues),
  sdBeta = apply(Bvalues, 1, sd),
  meanM  = rowMeans(Mvalues),
  sdM = apply(Mvalues, 1, sd)
)

# Make plot
p <- ggplot() +
  geom_point(data = plot_all, aes(x = meanBeta, y = sdBeta), 
             color = "lightgrey", alpha = 0.5) +
  geom_point(data = plotDF, aes(x = meanBeta, y = sdBeta), 
             color = RColorBrewer::brewer.pal(3, "Dark2")[1], alpha = 0.5) +
  xlim(c(0,1)) +
  ylim(c(0,0.5)) +
  xlab("Mean of \u03b2-values") +
  ylab("Standard Deviation of \u03b2-values") +
  ggtitle(subtitle) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "none",
        legend.title = element_text())

# Save plot
ggsave(p, file = paste0("MeanVsSD_", FeatureSelection, ".png"), width = 8, height = 6)


# Get beta-values mean and SD of specific probe annotation
load("E:/Thesis/MLData/probe_annotation.RData")

plot_all<- inner_join(plot_all, probe_annotation, by = c("ID" = "ID"))

p_island <- ggplot() +
  geom_point(data = plot_all, aes(x = meanBeta, y = sdBeta, color = Relation_to_Island), alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  xlim(c(0,1)) +
  ylim(c(0,0.5)) +
  xlab("Mean of \u03b2-values") +
  ylab("Standard Deviation of \u03b2-values") +
  ggtitle("Relation to CpG Island") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "right",
        legend.title = element_blank())

ggsave(p_island, file = "MeanVsSD_island.png", width = 8, height = 6)


p_island <- ggplot() +
  geom_point(data = plot_all, aes(x = meanBeta, y = sdBeta, color = Class), alpha = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  xlim(c(0,1)) +
  ylim(c(0,0.5)) +
  xlab("Mean of \u03b2-values") +
  ylab("Standard Deviation of \u03b2-values") +
  ggtitle("Probe Location") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "right",
        legend.title = element_blank())

ggsave(p_island, file = "MeanVsSD_class.png", width = 8, height = 6)



###############################################################################

# 2. Multiple regression

###############################################################################


features <- rownames(X_nonTest_varM)

formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                  paste0("0 + ", paste(colnames(Y_nonTest[,7:12]), collapse = " + ")))

# Make data matrix

# Scale data matrix
#dataMatrix_scaled <- t((X_nonTest_varM - rowMeans(X_nonTest_varM))/(apply(X_nonTest_varM,1,sd)))
dataMatrix <- cbind(t(X_nonTest_varM),Y_nonTest[,7:12])
model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))

fittedValues <- fitted(model)
residualValues <- residuals(model)

sse <- colSums(residualValues^2)
ssr <- colSums(fittedValues^2)

Rsquared <- ssr/(ssr + sse)

cellFeatures <- names(Rsquared)[Rsquared > 0.9]

X_nonTest_cell <- X_nonTest_varM[setdiff(cellFeatures, probe_annotation$ID[probe_annotation$Chr == "chrX"]),]

p_cellType <- ggplot() +
  geom_point(data = plot_all, aes(x = meanBeta, y = sdBeta), 
             color = "lightgrey", alpha = 0.5) +
  geom_point(data = plot_all[plot_all$ID %in% cellFeatures,], aes(x = meanBeta, y = sdBeta, color ='"Cell Type"-probes'), 
             color = RColorBrewer::brewer.pal(6, "Dark2")[6], alpha = 0.5) +
  xlim(c(0,1)) +
  ylim(c(0,0.5)) +
  xlab("Mean of \u03b2-values") +
  ylab("Standard Deviation of \u03b2-values") +
  ggtitle('"Cell Type"-probes') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "none",
        legend.title = element_text())

ggsave(p_cellType, file = "MeanVsSD_cellType.png", width = 8, height = 6)



###############################################################################

# 3. Venn diagram

###############################################################################

library(ggvenn)

plotList <- list(rownames(X_nonTest_S), 
                 rownames(X_nonTest_var), 
                 rownames(X_nonTest_varM))

names(plotList) <- c("S-score", "Variance (\u03b2)", "Variance (M)")

p <- ggvenn(plotList, show_percentage = FALSE)
ggsave(p, file = "FeatureSelectionVenn.png", width = 6, height = 6)


###############################################################################

# 4. PCA

###############################################################################

# Perform PCA
pcaList <-  prcomp(t(Bvalues),        
                   retx = TRUE,
                   center =TRUE,
                   scale = TRUE,
                   rank. = 10)

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
PCAscores$Sex <- ifelse(PCAscores$Sex == 1, "Male", "Female")
p_12 <- ggplot(data = PCAscores_var, aes(x = PC1, y = PC2)) +
  stat_ellipse(geom = "polygon",
               aes(fill = factor(Sex)),
               type = "norm", 
               alpha = 0.25,
               level = 0.99) +
  geom_point(alpha = 0.9, size = 2, aes(color = Neu)) +
  xlab(paste0("PC1 (", explVar_var[1],"%)")) +
  ylab(paste0("PC2 (", explVar_var[2],"%)")) +
  ggtitle("Variance (\u03b2)-based Feature Selection") +
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
ggsave(p_12, file = paste0("PCAplot_1vs2_", FeatureSelection, ".png"), width = 8, height = 6)

#=============================================================================#
# PCA correlations
#=============================================================================#

# Get principal components
PCs <- PCAscores[,1:5]
colnames(PCs) <- paste0(colnames(PCs), " (", explVar[1:5], "%)")

# Get meta data
meta <- PCAscores[,c("Age", "Sex", colnames(cellType))]
meta <- meta[,-c(9,10)]
meta$Sex <- ifelse(meta$Sex == "Male", 0,1)
colnames(meta) <- c("Age", "Sex", "CD8 T-cells", "CD4 T-cells", "NK cells", "B-cells", "Monocytes",
                    "Neutrophils")

# Get correlations
corDF <- expand.grid(colnames(PCs), colnames(meta))
colnames(corDF) <- c("PCs", "Meta")
corDF$Correlation <- rep(NA, nrow(corDF))
for (i in 1:nrow(corDF)){
  corDF$Correlation[i] <- cor(PCs[,corDF$PCs[i]], meta[,corDF$Meta[i]], method = "spearman")
}

# Make plot
p <- ggplot() +
  geom_point(data = corDF, aes(x = PCs, y = Meta, color = Correlation, size = abs(Correlation))) +
  labs(color = "Spearman\nCorrelation", size = "|Spearman\nCorrelation|") +
  guides(size = "none") +
  ggtitle(subtitle) +
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

# save plot
ggsave(p,file = paste0("correlationPlot_", FeatureSelection, ".png"), width = 8, height = 6)


