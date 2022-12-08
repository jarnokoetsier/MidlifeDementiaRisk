# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load cell type composition and phenotype data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Load model information
Score = "CAIDE1"
FeatureSelection = "Non"
load(paste0("CV_", Score, "_", FeatureSelection,".RData"))

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

################################################################################

# Alpha vs Lambda

################################################################################

#*****************************************************************************#
# Heatmap
#*****************************************************************************#

# Make plot
RMSE_heatmap <- ggplot() +
  geom_tile(data = trainResults, 
            aes(x = log(lambda), y = alpha, fill = RMSE)) +
  geom_point(data = trainResults[(trainResults$alpha == optAlpha) & (trainResults$lambda == optLambda),], 
             aes(x = log(lambda), y = alpha), 
             color = 'red', size = 3) +
  geom_label_repel(data = trainResults[(trainResults$alpha == optAlpha) & (trainResults$lambda == optLambda),], 
                   aes(x = log(lambda), y = alpha, label = paste0("Alpha: ", round(alpha,2), "\n",
                                                                  "Lambda: ", round(lambda,2), "\n",
                                                                  "Mean RMSE: ", round(RMSE,2), "\n",
                                                                  "SD RMSE: ", round(RMSESD,2))), 
                   color = "black", fill = "white", alpha = 1) +
  xlab("log \u03bb") +
  ylab("\u03b1") +
  labs(fill = "Mean\nRMSE") +
  scale_y_continuous(breaks = alphaCV[c(TRUE, FALSE)]) +
  ggtitle(Score, subtitle = subtitle) +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(legend.title = element_text(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(RMSE_heatmap, file = paste0(Score, "_", FeatureSelection, "_RMSE_heatmap.png"), width = 8, height = 6)


################################################################################

# Observed vs Predicted

################################################################################

#*****************************************************************************#
# Scatter plot
#*****************************************************************************#

# Make plot
ObsVsPred <- ggplot(ObsPred_CV) +
  geom_bin2d(aes(x = Observed, y = Predicted), bins = 100) +
  ggtitle(Score, subtitle = subtitle) +
  geom_abline(aes(intercept = 0, slope = 1), color = "red", linetype = "dashed", linewidth = 1.5) +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

# Save plot
ggsave(ObsVsPred, file = paste0(Score, "_", FeatureSelection, "_ObsVsPred.png"), width = 8, height = 6)

# Get R-squared
R2(ObsPred_CV$Predicted, ObsPred_CV$Observed)

#*****************************************************************************#
# Box plot
#*****************************************************************************#

ObsPred_CV$Class <- rep("Intermediate Risk (4-7)",nrow(ObsPred_CV))
ObsPred_CV$Class[ObsPred_CV$Observed <= 3] <- "Low Risk (0-3)"
ObsPred_CV$Class[ObsPred_CV$Observed >= 8] <- "High Risk (8-14)"
ObsPred_CV$Class <- factor(ObsPred_CV$Class, levels = c("Low Risk (0-3)", 
                                                        "Intermediate Risk (4-7)",
                                                        "High Risk (8-14)"))

ObsVsPred_box <- ggplot(ObsPred_CV) +
  geom_jitter(aes(x = Class, y = Predicted, color = Class), alpha = 0.3) +
  geom_boxplot(aes(x = Class, y = Predicted, fill = Class), alpha = 0.5) +
  xlab(NULL) +
  ylab(paste0("Predicted", Score, "Score")) +
  ggtitle(Score, subtitle = subtitle) +
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Reds") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(ObsVsPred_box, file = paste0(Score, "_", FeatureSelection, "_ObsVsPred_box.png"), width = 8, height = 6)


################################################################################

# Regression Coefficients

################################################################################


#*****************************************************************************#
# Heatmap
#*****************************************************************************#

# Evaluate stability of regression coefficients
coefs_finalModel <- as.data.frame(as.matrix(coef(finalModel)))

# only look at features in final model
rownames(coefs) <- rownames(coefs_finalModel)
coefs1 <- coefs[coefs_finalModel[,1] != 0,]

# remove intercept
coefs1 <- coefs1[-1,]

# Format data for plotting
coefPlot <- gather(as.data.frame(coefs1))
coefPlot$cpg <- rep(rownames(coefs1), ncol(coefs1))
coefPlot$avgValue <- rep(rowMeans(coefs1), ncol(coefs1))
coefPlot$finalModel <- rep(coefs_finalModel[-c(1,which(coefs_finalModel[,1] == 0)),1], ncol(coefs1))

# Main plot: heatmap
main <- ggplot(coefPlot, aes(x = fct_reorder(cpg, avgValue), y = key, fill = value)) +
  geom_tile()+
  xlab("CpG sites") +
  ylab("Folds in\nrepeated CV") +
  labs(fill = "Regression \nCoefficient") +
  scale_fill_viridis_c(limits = c(-0.8, 0.8), oob = scales::squish) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom")

# Top plot: scatter
top <- ggplot(coefPlot) +
  #geom_point(aes(x = fct_reorder(cpg, avgValue), y = avgValue), color = "blue") +
  geom_point(aes(x = fct_reorder(cpg, avgValue), y = finalModel, color = finalModel)) +
  ylab("Coefficients\nFinal Model") +
  scale_color_viridis_c(limits = c(-0.8, 0.8), oob = scales::squish) +
  scale_fill_viridis_c(limits = c(-0.8, 0.8), oob = scales::squish) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# Combine plots
p <- ggpubr::ggarrange(top,
                       main,
                       heights = c(2,8), 
                       nrow = 2,
                       ncol = 1,
                       align = "v")

ggsave(p, file = paste0(Score, "_", FeatureSelection, "_coefPlot.png"), width = 7, height = 8)



#*****************************************************************************#
# PCA on selected features
#*****************************************************************************#

# Load data matrix of selected features
load("X_coefs.RData")

# Perform PCA
pcaList <-  prcomp(t(X_coefs),        
                   retx = TRUE,
                   center =TRUE,
                   scale = TRUE,
                   rank. = 10)

# Get PCA scores
PCAscores <- as.data.frame(pcaList$x)
PCAscores$ID <- rownames(PCAscores)

# Combine with sample info
PCAscores <- inner_join(PCAscores, dat, by = c("ID" = "Basename"))

# Combine with cell type composition
PCAscores <- inner_join(PCAscores,cellType, by = c("ID" = "ID"))

# Combine with scores
if (Score == "CAIDE1") {
  load("E:/Thesis/EXTEND/Phenotypes/CAIDE.Rdata")
  PCAscores <- inner_join(PCAscores,CAIDE, by = c("ID" = "Basename"))
}

if (Score == "LIBRA") {
  load("E:/Thesis/EXTEND/Phenotypes/EPILIBRA.Rdata")
  PCAscores <- inner_join(PCAscores,EPILIBRA, by = c("ID" = "Basename"))
}

# Get explained variance
explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)


# Plot PC1 vs PC2
p_12 <- ggplot(data = PCAscores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.9, size = 2, aes(color = CAIDE)) +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  labs(color = "CAIDE1") +
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
ggsave(p,file = paste0(Score, "_", FeatureSelection, "_PC1vPC2.png"), width = 8, height = 6)


# Get Principal Components
PCs <- PCAscores[,1:5]
colnames(PCs) <- paste0(colnames(PCs), " (", explVar[1:5], "%)")

# Get cell type composition
meta_cell <- PCAscores[,colnames(cellType)[-c(7,8)]]
colnames(meta_cell) <- c("CD8 T-cells", "CD4 T-cells", "NK cells", "B-cells", "Monocytes",
                         "Neutrophils")

# Get CAIDE factors
meta_caide <- PCAscores[,176:183]
colnames(meta_caide) <- c("Age", "Sex", "Education", "Blood Pressure", "BMI",
                          "Serum HDL", "Physical Activity", "CAIDE1")

# Combine cell type composition and CAIDE factors
meta <- cbind.data.frame(meta_caide, meta_cell)

# Get correlations
corDF <- expand.grid(colnames(PCs), colnames(meta))
colnames(corDF) <- c("PCs", "Meta")
corDF$Correlation <- rep(NA, nrow(corDF))
for (i in 1:nrow(corDF)){
  corDF$Correlation[i] <- cor(PCs[,corDF$PCs[i]], meta[,corDF$Meta[i]], method = "spearman")
}

# Separate cell type correlations from CAIDE1 factors correlations
corDF$Class <- ifelse(corDF$Meta %in% colnames(meta_cell), "Cell Types", "CAIDE1 Factors")

# Make plot
p <- ggplot() +
  geom_point(data = corDF, aes(x = PCs, y = Meta, color = Correlation, size = abs(Correlation))) +
  facet_grid(rows =  vars(Class), scales = "free", space = "free") +
  labs(color = "Spearman\nCorrelation", size = "|Spearman\nCorrelation|") +
  guides(size = "none") +
  ggtitle("Features in Final Model") +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        strip.background = element_rect(
          color="black", fill="#1B9E77", linewidth = 1.5, linetype="solid"
        ))

# Save plot
ggsave(p,file = "correlationPlot_coefsModel.png", width = 8, height = 6)


#*****************************************************************************#
# Multiple regression
#*****************************************************************************#


#=============================================================================#
# Cell Types
#=============================================================================#

# Multiple regression
X_CAIDE1_coefs <- X_coefs[,Y_CAIDE1$Basename] 
features <- rownames(X_CAIDE1_coefs)
cellType_coefs <- cellType[colnames(X_CAIDE1_coefs),]
formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                  paste0("0 + ", paste(colnames(cellType_coefs[,1:6]), collapse = " + ")))

# Add factors
X_coefs_scaled <- (X_CAIDE1_coefs - rowMeans(X_CAIDE1_coefs))/(apply(X_CAIDE1_coefs, 1,sd))
dataMatrix <- cbind(t(X_coefs_scaled),cellType_coefs[,1:6])

# Fit model
model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
coeff <- coef(model)

fittedValues <- fitted(model)
residualValues <- residuals(model)

sse <- colSums(residualValues^2)
ssr <- colSums(fittedValues^2)

plotDF <- data.frame(Rsquared = ssr/(ssr + sse),
                     Probe = names(ssr))

ggplot(plotDF) +
  geom_histogram(aes(x = Rsquared), color = "white", bins = 20) +
  xlab(expression(R^2))+
  ylab("Count") +
  theme_classic()


#=============================================================================#
# CAIDE factors
#=============================================================================#

# Multiple regression 
X_CAIDE1_coefs <- X_coefs[,Y_CAIDE1$Basename] 
features <- rownames(X_CAIDE1_coefs)
formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                  paste0("0 + ", paste(colnames(Y_CAIDE1[,14:20]), collapse = " + ")))

# Add factors
X_coefs_scaled <- (X_CAIDE1_coefs - rowMeans(X_CAIDE1_coefs))/(apply(X_CAIDE1_coefs, 1,sd))
dataMatrix <- cbind(t(X_coefs_scaled),Y_CAIDE1[,14:20])

# Fit model
model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))

fittedValues <- fitted(model)
residualValues <- residuals(model)

sse <- colSums(residualValues^2)
ssr <- colSums(fittedValues^2)

plotDF <- data.frame(Rsquared = ssr/(ssr + sse),
                     Probe = names(ssr))

ggplot(plotDF) +
  geom_histogram(aes(x = Rsquared), color = "white", bins = 20) +
  xlab(expression(R^2))+
  ylab("Count") +
  theme_classic()

###############################################################################

# Mean Beta-value vs SD Beta-value

###############################################################################

# Load all beta/M and their sd
load("plot_all.RData")

# Load betas of probes in final model
load("X_coefs.RData")

# Format data
plotCoefs <- data.frame(
  ID = rownames(X_coefs),
  meanBeta = rowMeans(X_coefs),
  sdBeta = apply(X_coefs, 1, sd),
  meanM  = rowMeans(Mvalues_coefs),
  sdM = apply(Mvalues_coefs, 1, sd)
)

plot_all$ID <- rownames(plot_all)


# Make plot
p_coefs <- ggplot() +
  geom_point(data = plot_all, aes(x = meanBeta, y = sdBeta), 
             color = "lightgrey", alpha = 0.5) +
  geom_point(data = plotCoefs, aes(x = meanBeta, y = sdBeta), 
             color = RColorBrewer::brewer.pal(4, "Dark2")[4], alpha = 0.5) +
  xlim(c(0,1)) +
  ylim(c(0,0.5)) +
  xlab("Mean of \u03b2-values") +
  ylab("Standard Deviation of \u03b2-values") +
  ggtitle("Features in Final Model") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "none",
        legend.title = element_text())

# Save plot
ggsave(p_coefs, file = "MeanVsSD_coefs.png", width = 8, height = 6)




################################################################################

# Compare performances of feature selection approaches

################################################################################


# Score and feature selection method
Score = "CAIDE1"
methods = c("S", "var", "varM", "Non", "KS")

# Retrieve performance in CV for the different feature selection methods
Performances <- list()
for (i in 1:length(methods)){
  
  # load output
  load(paste0("CV_", Score, "_", methods[i],".RData"))
  
  # Put performance into list
  Performances[[i]] <- perf
}
names(Performances) <- methods

# Combine into single data frame
performanceDF <- data.frame(
  RMSE = unlist(Performances),
  FeatureSelection = c(rep("S-score-based Feature Selection", 25),
                       rep("Variance (\u03b2)-based Feature Selection", 25),
                       rep("Variance (M)-based Feature Selection", 25),
                       rep("No Feature Selection", 25),
                       rep("Kennard-Stone-like Feature Selection", 25)),
)

# Make plot
p <- ggplot(performanceDF) +
  geom_boxplot(aes(x = FeatureSelection, y = RMSE, fill = FeatureSelection), alpha = 0.3) +
  geom_point(aes(x = FeatureSelection, y = RMSE, color = FeatureSelection), 
             position=position_jitterdodge(jitter.width = 0.2), size = 2) +
  xlab("Feature Selection Method") +
  ggtitle(Score) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(ObsVsPred, file = paste0(Score, "_RMSE_boxplot.png"), width = 8, height = 6)

