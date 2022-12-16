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
files <- list.files(paste0("X_", FeatureSelection))
for (f in files){
  load(paste0("X_", FeatureSelection, "/",f))
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
if (FeatureSelection == "varCor"){
  subtitle = "Variance (\u03b2, Cor)-based Feature Selection"
}
if (FeatureSelection == "varMCor"){
  subtitle = "Variance (M, Cor)-based Feature Selection"
}

###############################################################################

# 1. Mean Beta-value vs SD Beta-value

###############################################################################

# Load all beta/M and their sd
load("plot_all.RData")

# Get beta- and M-values
Bvalues <- X_nonTest_varCor
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
             color = RColorBrewer::brewer.pal(5, "Dark2")[5], alpha = 0.5) +
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




methods <- c("S-score", "Variance (\u03b2)", "Variance (M)", "Variance (\u03b2, Cor)",
             "Variance (M, Cor)")
plotDF = NULL
for (i in 1:length(methods)){
  # Get data
  if (methods[i] == "S-score"){
    dataMatrix <- X_nonTest_S
  }
  if (methods[i] == "Variance (\u03b2)"){
    dataMatrix <- X_nonTest_var
  }
  if (methods[i] == "Variance (M)"){
    dataMatrix <- X_nonTest_varM
  }
  if (methods[i] == "Variance (\u03b2, Cor)"){
    dataMatrix <- X_nonTest_varCor
  }
  
  if (methods[i] == "Variance (M, Cor)"){
    dataMatrix <- X_nonTest_varMCor
  }
  
  
  # Get Features
  features <- rownames(dataMatrix)
  
  # Make formula
  formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                    paste0("0 + ", paste(colnames(Y_nonTest[,7:12]), collapse = " + ")))
  
  # Scale data matrix
  dataMatrix_scaled <- t((dataMatrix - rowMeans(dataMatrix))/(apply(dataMatrix,1,sd)))
  
  # Combine with cell type composition
  dataMatrix <- cbind(dataMatrix_scaled,Y_nonTest[,7:12])
  
  # Make linear model
  model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
  
  # Get fitted values
  fittedValues <- fitted(model)
  
  # Get residual values
  residualValues <- residuals(model)
  
  # Calculate R-squared
  sse <- colSums(residualValues^2)
  ssr <- colSums(fittedValues^2)
  Rsquared <- ssr/(ssr + sse)
  
  # Prepare data for plotting
  temp <- data.frame(
    Value = c(sum(Rsquared > 0.1),
              sum(Rsquared > 0.5),
              sum(Rsquared > 0.8)),
    FeatureSelection = c(rep(methods[i], 3)),
    Threshold = rep(c("R2 > 0.1","R2 > 0.5","R2 > 0.8"))
  )
  
  plotDF <- rbind.data.frame(plotDF, temp)
}
plotDF$FeatureSelection <- factor(plotDF$FeatureSelection, levels = methods)

# Make plot
p <- ggplot() +
  geom_bar(data = plotDF, aes(x = FeatureSelection, y = Value, fill = Threshold),
           stat="identity", position=position_dodge(), color = "black", alpha = 0.8) +
  xlab("Feature Selection Method") +
  ylab("# Features") +
  scale_fill_manual(breaks = c("R2 > 0.1","R2 > 0.5","R2 > 0.8"),
                      labels = c(expression(R^2 > 0.1), 
                                 expression(R^2 > 0.5),
                                 expression(R^2 > 0.8)),
                      values = RColorBrewer::brewer.pal(3, "Reds")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")


# Save plot
ggsave(p, file = "R2_cellType.png", width = 8, height = 6)




###############################################################################

# 3. Venn diagram

###############################################################################

library(ggvenn)

plotList <- list(rownames(X_nonTest_S), 
                 rownames(X_nonTest_var), 
                 rownames(X_nonTest_varM))

names(plotList) <- c("S-score", "Variance (\u03b2)", "Variance (M)")

p <- ggvenn(plotList, show_percentage = FALSE, fill_color = RColorBrewer::brewer.pal(3, "Dark2"))
ggsave(p, file = "FeatureSelectionVenn.png", width = 6, height = 6)


###############################################################################

# 4. Heatmap of shared features

###############################################################################

featureList <- list(features_S = rownames(X_nonTest_S),
                    features_var = rownames(X_nonTest_var),
                    features_varM = rownames(X_nonTest_varM),
                    features_varCor = rownames(X_nonTest_varCor),
                    features_varMCor = rownames(X_nonTest_varMCor))
# Load all files
m <- c("S", "var", "varM", "varCor", "varMCor")

for (FeatureSelection in m){
  files <- list.files(paste0("X_", FeatureSelection))
  for (f in files){
    load(paste0("X_", FeatureSelection, "/",f))
  }
}


# Format data
combinations <- expand.grid(c(1:5),c(1:5),c(1:5),c(1:5),c(1:5))
colnames(combinations) <- m

for (i in 1:length(methods)){
  combinations[combinations == i] <- m[i]
}

combinations1 <- combinations
for (i in 1:ncol(combinations)) {
  combinations1[,i] <- ifelse(apply(combinations, 1, function(x) {m[i] %in% x}),1,0)
}
combinations1 <- unique(combinations1)


combinations1$Sum <- apply(combinations1, 1, function(x) {
  length(Reduce(intersect, featureList[which(x != 0)]))
})

rownames(combinations1) <- as.character(1:nrow(combinations1))
combinations1$ID <- as.character(1:nrow(combinations1))
colnames(combinations1) <- c("S-score", "Variance (\u03b2)", "Variance (M)",
                             "Variance (\u03b2, Cor)", "Variance (M, Cor)",
                             "Sum", "ID")

levelsID <- rev(names(sort(rowSums(combinations1[,1:5]))))
levelsSelection <- c("S-score", "Variance (\u03b2)", "Variance (M)",
                     "Variance (\u03b2, Cor)", "Variance (M, Cor)")


# Make heatmap part of the plot
plotDF <- gather(combinations1[,1:5])
plotDF$ID <- factor(rep(combinations1$ID,5),
                    levels = levelsID)
plotDF$key <- factor(plotDF$key,
                     levels = levelsSelection)

plotDF <- plotDF[plotDF$value == 1,]

main <- ggplot(plotDF) +
  geom_tile(aes(x = ID, y = key, fill = key),
            color = 'black', width = 0.9, height = 0.9, size = 0.7) +
  ylab("Feature Selection Method") +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 10, color = "black"),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())



# Make bar chart part of the plot
plotDF1 <- combinations1[,c("Sum", "ID")]
plotDF1$ID <- factor(plotDF1$ID,
                     levels = levelsID)

top <- ggplot(plotDF1) +
  geom_bar(aes(x = ID, y = Sum), stat = "identity", color = "black", fill = "grey", size = 0.7) +
  ylab("# Features") +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, size = 10, color = "black"),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90, face = "italic"))


library(patchwork)

p <- top + main +
  plot_layout(nrow = 2, ncol = 1) +
  plot_layout(heights = c(1,1))

ggsave(p, file = "UpsetDiagram.png", width = 10, height = 6)


features_S <- rownames(X_nonTest_S)
features_var <- rownames(X_nonTest_var)
features_varM <- rownames(X_nonTest_varM)
features_varCor <- rownames(X_nonTest_varCor)
features_varMCor <- rownames(X_nonTest_varMCor)

plotDF <- data.frame(
  Source = rep(c("S-score", "Variance (\u03b2)", "Variance (M)",
                 "Variance (\u03b2, Cor)", "Variance (M, Cor)"),5),
  Target =  c(rep("S-score",5), rep("Variance (\u03b2)",5), rep("Variance (M)",5),
              rep("Variance (\u03b2, Cor)",5), rep("Variance (M, Cor)",5)),
  value = c(length(intersect(features_S, features_S)),
            length(intersect(features_var, features_S)),
            length(intersect(features_varM, features_S)),
            length(intersect(features_varCor, features_S)),
            length(intersect(features_varMCor, features_S)),
            length(intersect(features_S, features_var)),
            length(intersect(features_var, features_var)),
            length(intersect(features_varM, features_var)),
            length(intersect(features_varCor, features_var)),
            length(intersect(features_varMCor, features_var)),
            length(intersect(features_S, features_varM)),
            length(intersect(features_var, features_varM)),
            length(intersect(features_varM, features_varM)),
            length(intersect(features_varCor, features_varM)),
            length(intersect(features_varMCor, features_varM)),
            length(intersect(features_S, features_varCor)),
            length(intersect(features_var, features_varCor)),
            length(intersect(features_varM, features_varCor)),
            length(intersect(features_varCor, features_varCor)),
            length(intersect(features_varMCor, features_varCor)),
            length(intersect(features_S, features_varMCor)),
            length(intersect(features_var, features_varMCor)),
            length(intersect(features_varM, features_varMCor)),
            length(intersect(features_varCor, features_varMCor)),
            length(intersect(features_varMCor, features_varMCor))
            )
)
levels = unique(plotDF$Source)


plotDF$Combined <- paste(plotDF$Source, plotDF$Target, sep = "_")
plotDF$value[plotDF$Combined %in% c(paste0(levels[2], "_", levels[1]),
                                    paste0(levels[3], "_", levels[1]),
                                    paste0(levels[4], "_", levels[1]),
                                    paste0(levels[5], "_", levels[1]),
                                    paste0(levels[3], "_", levels[2]),
                                    paste0(levels[4], "_", levels[2]),
                                    paste0(levels[5], "_", levels[2]),
                                    paste0(levels[4], "_", levels[3]),
                                    paste0(levels[5], "_", levels[3]),
                                    paste0(levels[5], "_", levels[4])
                                    )] <- NA

plotDF$Source <- factor(plotDF$Source, levels = levels)
plotDF$Target <- factor(plotDF$Target, levels = rev(levels))

p <- ggplot(plotDF) +
  geom_tile(aes(x = Source, y = Target, fill = value, color = ifelse(is.na(value), NA, "1")),
            width = 0.9, height = 0.9, alpha = 0.8) +
  geom_label(aes(x = Source, y = Target, label = value), 
             alpha = 0.5, label.size = NA) +
  scale_fill_viridis_c(na.value = "white") +
  scale_color_manual(values = "black", na.value = "white") +
  theme_void() +
  theme(legend.position = "none",
        axis.text.y = element_text(hjust = 1),
        axis.text.x = element_text())

ggsave(p, file = "FeatureSelectionHeatmap.png", width = 8, height = 6)

###############################################################################

# 4. PCA explained variances

###############################################################################

methods <- c("S-score", "Variance (\u03b2)", "Variance (M)", "Variance (\u03b2, Cor)",
             "Variance (M, Cor)")
plotExplVar = NULL
for (i in 1:length(methods)){
  # Get data
  if (methods[i] == "S-score"){
    dataMatrix <- X_nonTest_S
  }
  if (methods[i] == "Variance (\u03b2)"){
    dataMatrix <- X_nonTest_var
  }
  if (methods[i] == "Variance (M)"){
    dataMatrix <- X_nonTest_varM
  }
  if (methods[i] == "Variance (\u03b2, Cor)"){
    dataMatrix <- X_nonTest_varCor
  }
  
  if (methods[i] == "Variance (M, Cor)"){
    dataMatrix <- X_nonTest_varMCor
  }
  
  # Perform PCA
  pcaList <-  prcomp(t(dataMatrix),        
                     retx = TRUE,
                     center =TRUE,
                     scale = TRUE)
  
  # Get explained variance
  explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)
  cumVar <- sapply(1:length(explVar),function(x){sum(explVar[1:x])})

  temp <- data.frame(
    explVar = explVar,
    cumVar = cumVar,
    FeatureSelection = rep(methods[i], length(explVar))
  )
  
  plotExplVar <- rbind.data.frame(plotExplVar, temp)
}
plotExplVar$PC <- factor(rep(paste0("PC", 1:924), length(methods)),
                             levels = paste0("PC", 1:924))

plotExplVar$FeatureSelection <- factor(plotExplVar$FeatureSelection,
                                       levels = methods)

p  <- ggplot(plotExplVar) +
  geom_step(aes(x = PC, y = cumVar, group = FeatureSelection, color = FeatureSelection),
            linewidth = 1.5) +
  xlab("Principal Components (PC1 - 924)") +
  ylab("Cumulative Explained Variance (%)") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.title  = element_blank(),
        axis.text.x = element_blank())#element_text(angle = 45, vjust = 0.5))
  
ggsave(p, file = "FeatureSelectionExplVar.png", width = 8, height = 5)

###############################################################################

# 5. PCA

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
p_12 <- ggplot(data = PCAscores, aes(x = PC1, y = PC2)) +
  stat_ellipse(geom = "polygon",
               aes(fill = factor(Sex)),
               type = "norm", 
               alpha = 0.25,
               level = 0.99) +
  geom_point(alpha = 0.9, size = 2, aes(color = Neu)) +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
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


