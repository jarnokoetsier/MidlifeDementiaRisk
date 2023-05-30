# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(ggbreak)

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

###############################################################################

# 1. Mean Beta-value vs SD Beta-value

###############################################################################

# Feature selection data
FeatureSelection = "Cor_CAIDE1"
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/",f))
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
if (FeatureSelection == "PC"){
  subtitle = "PCA-based Feature Selection"
}
if (FeatureSelection == "KS"){
  subtitle = "Kennard-Stone-like Feature Selection"
}
if (FeatureSelection == "Cor_CAIDE1"){
  subtitle = "Correlation-based Feature Selection"
}
if (FeatureSelection == "CorKS"){
  subtitle = "Correlation + KS-based Feature Selection"
}

# Load all beta/M and their sd
load("plot_all.RData")

# Get beta- and M-values
Bvalues <- X_nonTest_Cor
Mvalues <- log2(Bvalues/(1 - Bvalues))

# Format data
plotDF <- data.frame(
  ID = rownames(Bvalues),
  meanB = rowMeans(Bvalues),
  sdB = apply(Bvalues, 1, sd),
  meanM  = rowMeans(Mvalues),
  sdM = apply(Mvalues, 1, sd)
)

# Make plot
p <- ggplot() +
  geom_point(data = plot_all, aes(x = meanB, y = sdB), 
             color = "lightgrey", alpha = 0.5) +
  geom_point(data = plotDF, aes(x = meanB, y = sdB), 
             color = RColorBrewer::brewer.pal(8, "Dark2")[7], alpha = 0.5) +
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



###############################################################################

# 2. % explained by cell type composition

###############################################################################


# Load all files
m <- c("S", "var", "varM", "varCor", "varMCor","KS", "Cor_CAIDE1")

for (FeatureSelection in m){
  files <- list.files(paste0("X/X_", FeatureSelection))
  for (f in files){
    load(paste0("X/X_", FeatureSelection, "/",f))
  }
}
methods <- c("S-score", "Variance (\u03b2)", "Variance (M)", "Variance (\u03b2, Cor)",
             "Variance (M, Cor)","KS-like", "Correlation")


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
  if (methods[i] == "PCA"){
    dataMatrix <- t(X_nonTest_PC)
  }
  if (methods[i] == "KS-like"){
    dataMatrix <- X_nonTest_KS
  }
  if (methods[i] == "Correlation"){
    dataMatrix <- X_nonTest_Cor
  }
  
  # Get Features
  features <- rownames(dataMatrix)
  
  # Make formula
  formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                    paste0("0 + ", paste(colnames(Y_nonTest[,7:12]), collapse = " + ")))
  
  # Scale data matrix
  dataMatrix_scaled <- t((dataMatrix - rowMeans(dataMatrix)))
  
  # Combine with cell type composition
  dataMatrix1 <- cbind(dataMatrix_scaled,Y_nonTest[,7:12])
  
  # Make linear model
  model <- lm(as.formula(formula), data = as.data.frame(dataMatrix1))
  
  # Get fitted values
  fittedValues <- fitted(model)
  
  # Get residual values
  residualValues <- residuals(model)
  
  # Calculate R-squared
  sse <- colSums(residualValues^2)
  ssr <- colSums(fittedValues^2)
  sst <- colSums(dataMatrix_scaled^2)
  explVar <- sum(ssr)/sum(ssr + sse)
  
  # Prepare data for plotting
  temp <- data.frame(
    Value = explVar,
    FeatureSelection = methods[i]
  )
  
  plotDF <- rbind.data.frame(plotDF, temp)
}
plotDF$FeatureSelection <- factor(plotDF$FeatureSelection, levels = methods)

# Make plot
p <- ggplot() +
  geom_bar(data = plotDF, aes(x = FeatureSelection, y = Value*100, fill = FeatureSelection),
           stat="identity", position=position_dodge(), color = "black", alpha = 0.8) +
  coord_flip() +
  xlab("Feature Selection Method") +
  ylab("Variance explained (%) by cell type composition") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none")


# Save plot
ggsave(p, file = "R2_cellType_all.png", width = 8, height = 6)


###############################################################################

# 2. Multiple regression

###############################################################################


# Load all files
m <- c("S", "var", "varM", "varCor", "varMCor","KS", "Cor_CAIDE1", "PC")

for (FeatureSelection in m){
  files <- list.files(paste0("X/X_", FeatureSelection))
  for (f in files){
    load(paste0("X/X_", FeatureSelection, "/",f))
  }
}

methods <- c("S-score", "Variance (\u03b2)", "Variance (M)", "Variance (\u03b2, Cor)",
             "Variance (M, Cor)", "KS-like", "Correlation", "PCA")

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
  if (methods[i] == "PCA"){
    dataMatrix <- t(X_nonTest_PC)
  }
  if (methods[i] == "KS-like"){
    dataMatrix <- X_nonTest_KS
  }
  if (methods[i] == "Correlation"){
    dataMatrix <- X_nonTest_Cor
  }
  
  # Get Features
  features <- rownames(dataMatrix)
  
  # Make formula
  formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                    paste0("0 + ", paste(colnames(Y_nonTest[,7:12]), collapse = " + ")))
  
  # Scale data matrix
  dataMatrix_scaled <- t((dataMatrix - rowMeans(dataMatrix))/(apply(dataMatrix,1,sd)))
  
  # Combine with cell type composition
  dataMatrix1 <- cbind(dataMatrix_scaled,Y_nonTest[,7:12])
  
  # Make linear model
  model <- lm(as.formula(formula), data = as.data.frame(dataMatrix1))
  
  # Get fitted values
  fittedValues <- fitted(model)
  
  # Get residual values
  residualValues <- residuals(model)
  
  # Calculate R-squared
  sse <- colSums(residualValues^2)
  ssr <- colSums(fittedValues^2)
  Rsquared <- ssr/(colSums(dataMatrix_scaled^2))
  
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
  geom_bar(data = plotDF[plotDF$FeatureSelection != "PCA",], aes(x = FeatureSelection, y = Value, fill = Threshold),
           stat="identity", position=position_dodge(), color = "black", alpha = 1) +
  #xlab("Feature Selection Method") +
  scale_y_break(c(1800,4500), ticklabels=c(4500,5000)) +
  xlab(NULL) +
  ylab("# Features") +
  coord_flip() +
  scale_fill_manual(breaks = c("R2 > 0.1","R2 > 0.5","R2 > 0.8"),
                      labels = c(expression(R^2 > 0.1), 
                                 expression(R^2 > 0.5),
                                 expression(R^2 > 0.8)),
                      values = RColorBrewer::brewer.pal(3, "Reds")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "right")


# Save plot
ggsave(p, file = "R2_cellType_horizontal.png", width = 8, height = 4)



###############################################################################

# 4. Upset diagram (overlap of features)

###############################################################################

# Load all files
m <- c("S", "var", "varM", "varCor", "varMCor","KS", "Cor_CAIDE1")

for (FeatureSelection in m){
  files <- list.files(paste0("X/X_", FeatureSelection))
  for (f in files){
    load(paste0("X/X_", FeatureSelection, "/",f))
  }
}

methods <- c("S-score", "Variance (\u03b2)", "Variance (M)", "Variance (\u03b2, Cor)",
             "Variance (M, Cor)", "KS-like", "Correlation")


featureList <- list(features_S = rownames(X_nonTest_S),
                    features_var = rownames(X_nonTest_var),
                    features_varM = rownames(X_nonTest_varM),
                    features_varCor = rownames(X_nonTest_varCor),
                    features_varMCor = rownames(X_nonTest_varMCor),
                    features_KS = rownames(X_nonTest_KS),
                    features_Cor = rownames(X_nonTest_Cor))

test <- data.frame(Features = unlist(featureList),
                   Method = rep(methods, each = 10000))

#install.packages("ggupset")
library(ggupset)
test <- test %>%
  as_tibble() %>%
  group_by(Features) %>%
  summarize(Method = list(Method))

p <- ggplot(data = test, aes(x=Method)) +
  geom_bar(fill = "grey", color = "black") +
  scale_x_upset() +
  xlab(NULL) +
  ylab("Count") +
  theme_classic() +
  theme_combmatrix()

ggsave(p, file = "upsetDiagram.png", width = 10, height = 6)



###############################################################################

# 4. PCA explained variances

###############################################################################

# Load all files
m <- c("S", "var", "varM", "varCor", "varMCor","KS", "Cor_CAIDE1", "PC")

for (FeatureSelection in m){
  files <- list.files(paste0("X/X_", FeatureSelection))
  for (f in files){
    load(paste0("X/X_", FeatureSelection, "/",f))
  }
}

methods <- c("S-score", "Variance (\u03b2)", "Variance (M)", "Variance (\u03b2, Cor)",
             "Variance (M, Cor)", "KS-like", "Correlation", "PCA")

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
  if (methods[i] == "KS-like"){
    dataMatrix <- X_nonTest_KS
  }
  if (methods[i] == "PCA"){
    dataMatrix <- t(X_nonTest_PC)
  }
  if (methods[i] == "Correlation"){
    dataMatrix <- X_nonTest_Cor
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

load("FeatureSelection/explVar_None.RData")
plotExplVar <- rbind.data.frame(plotExplVar, plotDF_PC)

plotExplVar$FeatureSelection <- factor(plotExplVar$FeatureSelection,
                                       levels = c(methods, "None"))

colors <- c(RColorBrewer::brewer.pal(n = 8, "Dark2"), "red")
p  <- ggplot(plotExplVar) +
  geom_step(aes(x = PC, y = cumVar, group = FeatureSelection, color = FeatureSelection),
            linewidth = 1.5) +
  geom_hline(yintercept = 100,  linetype = "dashed", linewidth = 1.5) +
  xlab("Principal Components (PC1 - 924)") +
  ylab("Cumulative Explained Variance (%)") +
  scale_color_manual(values = colors) +
  ylim(c(0,100)) +
  theme_classic() +
  theme(legend.title  = element_blank(),
        axis.text.x = element_blank())#element_text(angle = 45, vjust = 0.5))
  
ggsave(p, file = "FeatureSelectionExplVar.png", width = 10, height = 5)

# zoom-in
selected <- paste0("PC",1:10)
p  <- ggplot(plotExplVar[plotExplVar$PC %in% selected,]) +
  geom_step(aes(x = PC, y = cumVar, group = FeatureSelection, color = FeatureSelection),
            linewidth = 1.5) +
  xlab("Principal Components (PC1 - 10)") +
  ylab("Cumulative Explained Variance (%)") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.title  = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "right")#element_text(angle = 45, vjust = 0.5))

ggsave(p, file = "FeatureSelectionExplVar_zoom.png", width = 10, height = 5)


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


