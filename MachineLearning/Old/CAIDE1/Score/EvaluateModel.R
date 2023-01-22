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
FeatureSelection = "Cor"
load(paste0("CV_CAIDE1/CV_", Score, "_", FeatureSelection,"_EN.RData"))

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
if (FeatureSelection == "cor"){
  subtitle = "Correlation-based Feature Selection"
}
if (FeatureSelection == "KS"){
  subtitle = "Kennard-Stone-like Feature Selection"
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
  scale_y_continuous(breaks = seq(0.1,1,length.out = 10)[c(TRUE, FALSE)]) +
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
  ylab(paste0("Predicted ", Score, " Score")) +
  ggtitle(Score, subtitle = subtitle) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, "Reds")[2:4]) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 4, "Reds")[2:4]) +
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
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.8,0.8), oob = scales::squish) +
  #scale_fill_viridis_c(limits = c(-0.8, 0.8), oob = scales::squish) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "bottom")

# Top plot: scatter

# coefficients in final model
coefs_finalModel <- as.data.frame(as.matrix(coef(finalModel)))
names_coefs <- rownames(coefs_finalModel)[coefs_finalModel != 0][-1]
values_coefs <- coefs_finalModel[coefs_finalModel != 0,1][-1]
finalCoefs <- data.frame(CpG = names_coefs,
                         coefValue = values_coefs)
# Make plot
top <- ggplot(finalCoefs) +
  geom_bar(aes(x = fct_reorder(CpG, values_coefs), y = values_coefs, fill = values_coefs), stat = "identity", color = "black") +
  #geom_point(aes(x = fct_reorder(cpg, avgValue), y = finalModel, color = finalModel)) +
  ylab("Coefficients\nFinal Model") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.8,0.8), oob = scales::squish) +
  #scale_color_viridis_c(limits = c(-0.8, 0.8), oob = scales::squish) +
  #scale_fill_viridis_c(limits = c(-0.8, 0.8), oob = scales::squish) +
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
load("X_nonTest_coefs.RData")

# Perform PCA
pcaList <-  prcomp(t(X_nonTest_coeff),        
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
ggsave(p_12,file = paste0(Score, "_", FeatureSelection, "_PC1vPC2.png"), width = 8, height = 6)


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

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}
# coefficients in final model
coefs_finalModel <- as.data.frame(as.matrix(coef(finalModel)))
names_coefs <- rownames(coefs_finalModel)[coefs_finalModel != 0][-1]
values_coefs <- coefs_finalModel[coefs_finalModel != 0,1][-1]
finalCoefs <- data.frame(CpG = names_coefs,
                         coefValue = values_coefs)

# Check whether samples are in correct order
load("X/X_Cor/X_nonTest_Cor.RData")
X_nonTest_coeff <- X_nonTest_Cor[finalCoefs$CpG,]
all(Y_nonTest$Basename == colnames(X_nonTest_coeff))

# Get features
features <- rownames(X_nonTest_coeff)

# Make formula
formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                  paste0("0 + ", paste(colnames(Y_nonTest[,7:12]), collapse = " + ")))

# Mean center the data
X_coefs_scaled <- (X_nonTest_coeff - rowMeans(X_nonTest_coeff))

# Combine with dependent variables
dataMatrix <- cbind(t(X_coefs_scaled),Y_nonTest[,7:12])

# Fit model
model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))

# Get fitted and residual values
fittedValues <- fitted(model)
residualValues <- residuals(model)

# Calculate Rsquared
sse <- colSums(residualValues^2)
ssr <- colSums(fittedValues^2)
Rsquared = ssr/(ssr + sse)

# Calculate Cohen's f2 statistic
cohenF <- list()
cells <- colnames(Y_nonTest[,7:12])
for (i in 1:length(cells)){
  
  # Formula without cell type
  formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                    paste0("0 + ", paste(cells[-i], collapse = " + ")))
  
  # Fit model
  model_i <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
  
  # Get fitted and residual values
  fittedValues <- fitted(model_i)
  residualValues <- residuals(model_i)
  
  # Calculate Rsquared
  sse <- colSums(residualValues^2)
  ssr <- colSums(fittedValues^2)
  Rsquared_i <- ssr/(ssr + sse)
  
  # Calculate cohen's f2 statistic (local effect size)
  cohenF[[i]] <- (Rsquared - Rsquared_i)/(1-Rsquared)
}

# Get global effect size
globalEffect <- (Rsquared)/(1-Rsquared)

# Combine into data frame
cellNames <- c("CD8 T-cells", "CD4 T-cells", "NK cells", "B-cells", "Monocytes", "Neutrophils")
plotDF <- data.frame(cohenF = c(unlist(cohenF), globalEffect),
                     Effect = rep(c(cellNames, "Global"), each = nrow(X_coefs_scaled)),
                     CpG = rep(features,7))

# Reorder
plotDF$Effect <- factor(plotDF$Effect, levels = c(cellNames, "Global"))

# Combine with coefficient values in final model
plotDF <- inner_join(plotDF, finalCoefs, by = c("CpG" = "CpG"))
plotDF$CpG <- factor(plotDF$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))

# Make plot
main <- ggplot(plotDF) +
  geom_bar(aes(y = cohenF, x = CpG, fill = Effect), stat="identity") +
  facet_grid(rows = vars(Effect)) +
  xlab(NULL) +
  ylab(expression("Cohen's " ~ f^2)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        strip.background = element_rect(fill= "grey", linewidth = 0, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"))



finalCoefs$CpG <- factor(finalCoefs$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))
top <- ggplot(finalCoefs) +
  #geom_point(aes(x = fct_reorder(cpg, avgValue), y = avgValue), color = "blue") +
  geom_bar(aes(x = CpG, y = coefValue, fill = coefValue), stat = "identity", color = "black") +
  ylab("Coefficients\nFinal Model") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-0.8,0.8), oob = scales::squish) +
  scale_color_viridis_c(limits = c(-0.8, 0.8), oob = scales::squish) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

library(patchwork)
p <- top / main +
  plot_layout(heights = c(1,7))

# Save plot
ggsave(p,file = "cohenF_coefsModel.png", width = 8, height = 8)



# Distribution of R-squared
plotDF1 <- data.frame(Rsquared = Rsquared,
                     Probe = names(Rsquared))

p <- ggplot(plotDF1) +
  geom_histogram(aes(x = Rsquared), color = "white", bins = 20) +
  xlab(expression(R^2))+
  ylab("Count") +
  theme_classic()

ggsave(p,file = "R2_coefsModel.png", width = 8, height = 6)

#=============================================================================#
# CAIDE factors
#=============================================================================#

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}
# coefficients in final model
coefs_finalModel <- as.data.frame(as.matrix(coef(finalModel)))
names_coefs <- rownames(coefs_finalModel)[coefs_finalModel != 0][-1]
values_coefs <- coefs_finalModel[coefs_finalModel != 0,1][-1]
finalCoefs <- data.frame(CpG = names_coefs,
                         coefValue = values_coefs)

# Check whether samples are in correct order
X_CAIDE1_coeff <- X_nonTest_coeff[,Y_CAIDE1$Basename]

# Get features
features <- rownames(X_CAIDE1_coeff)

# Make formula
formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                  paste0("0 + ", paste(colnames(Y_CAIDE1[,14:20]), collapse = " + ")))

# Mean center the data
X_coefs_scaled <- (X_CAIDE1_coeff - rowMeans(X_CAIDE1_coeff))

# Combine with dependent variables
dataMatrix <- cbind(t(X_coefs_scaled),Y_CAIDE1[,14:20])

# Fit model
model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))

# Get fitted and residual values
fittedValues <- fitted(model)
residualValues <- residuals(model)

# Calculate Rsquared
sse <- colSums(residualValues^2)
ssr <- colSums(fittedValues^2)
Rsquared = ssr/(ssr + sse)

# Get global effect size
globalEffect <- (Rsquared)/(1-Rsquared)

# Calculate Cohen's f2 statistic
cohenF <- list()
factors <- colnames(Y_CAIDE1[,14:20])
for (i in 1:length(factors)){
  
  # Formula without factor
  formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                    paste0("0 + ", paste(factors[-i], collapse = " + ")))
  
  # Fit model
  model_i <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
  
  # Get fitted and residual values
  fittedValues <- fitted(model_i)
  residualValues <- residuals(model_i)
  
  # Calculate Rsquared
  sse <- colSums(residualValues^2)
  ssr <- colSums(fittedValues^2)
  Rsquared_i <- ssr/(ssr + sse)
  
  # Calculate cohen's f2 statistic (local effect size)
  cohenF[[i]] <- ((Rsquared - Rsquared_i)/(1-Rsquared))#/globalEffect
}



# Combine into data frame
factorNames <- c("Age", "Sex", "Edu", "BP", "BMI", "Chol", "Activity")
plotDF <- data.frame(cohenF = c(unlist(cohenF), globalEffect),
                     Effect = rep(c(factorNames, "Global"), each = nrow(X_coefs_scaled)),
                     CpG = rep(features,length(factorNames) +1))

# Reorder
plotDF$Effect <- factor(plotDF$Effect, levels = c(factorNames, "Global"))

# Combine with coefficient values in final model
plotDF <- inner_join(plotDF, finalCoefs, by = c("CpG" = "CpG"))
plotDF$CpG <- factor(plotDF$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))

# Make plot
main <- ggplot(plotDF) +
  geom_bar(aes(y = cohenF, x = CpG, fill = Effect), stat="identity") +
  facet_grid(rows = vars(Effect)) +
  xlab(NULL) +
  ylab(expression("Cohen's " ~ f^2)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        strip.background = element_rect(fill= "grey", linewidth = 0, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"))



finalCoefs$CpG <- factor(finalCoefs$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))
top <- ggplot(finalCoefs) +
  #geom_point(aes(x = fct_reorder(cpg, avgValue), y = avgValue), color = "blue") +
  geom_bar(aes(x = CpG, y = coefValue, fill = coefValue), stat = "identity", color = "black") +
  ylab("Coefficients\nFinal Model") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.8,0.8), oob = scales::squish) +
  scale_color_viridis_c(limits = c(-0.8, 0.8), oob = scales::squish) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

library(patchwork)
p <- top / main +
  plot_layout(heights = c(1,7))

# Save plot
ggsave(p,file = "cohenF_coefsModel_CAIDE1.png", width = 8, height = 8)



# Distribution of R-squared
plotDF1 <- data.frame(Rsquared = Rsquared,
                      Probe = names(Rsquared))

p <- ggplot(plotDF1) +
  geom_histogram(aes(x = Rsquared), color = "white", bins = 20) +
  xlab(expression(R^2))+
  ylab("Count") +
  theme_classic()

ggsave(p,file = "R2_coefsModel_CAIDE1.png", width = 8, height = 6)


###############################################################################

# Mean Beta-value vs SD Beta-value

###############################################################################

# Load all beta/M and their sd
load("plot_all.RData")

# Load betas of probes in final model
load("X_nonTest_coefs.RData")

Mvalues_coefs <- log2(X_nonTest_coeff/(1-X_nonTest_coeff))
# Format data
plotCoefs <- data.frame(
  ID = rownames(X_nonTest_coeff),
  meanBeta = rowMeans(X_nonTest_coeff),
  sdBeta = apply(X_nonTest_coeff, 1, sd),
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

# Compare performances of feature selection approaches (including test set)

################################################################################

# Score and feature selection method
Score = "CAIDE1"
methods = c("S", "var", "varM", "varCor", "varMCor","PC", "KS", "Non")

# Load data
for (i in 1:length(methods)){
    load(paste0("X/X_", methods[i], "/", "X_test_", methods[i], ".RData"))
}
# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}


# Retrieve performance in CV for the different feature selection methods
Performances <- list()
Performances_test <- list()
for (i in 1:length(methods)){
  
  # load output
  load(paste0("CV_CAIDE1/CV_", Score, "_", methods[i],".RData"))
  
  # Put performance into list
  Performances[[i]] <- perf
  
  # Performance in test data
  X_test <- get(paste0("X_test_", methods[i]))
  if (methods[i] != "PC"){
    testData <- log2(X_test/(1-X_test))
  }
  if (methods[i] == "PC"){
    testData <- t(X_test)
  }
  pred_test <- predict(finalModel, t(testData))
  Performances_test[[i]] <- RMSE(pred = pred_test,
                    obs = Y_test$CAIDE)
}
names(Performances) <- methods
names(Performances_test) <- methods



# Combine into single data frame
performanceDF <- data.frame(
  RMSE = unlist(Performances),
  FeatureSelection = c(rep("S-score", 25),
                       rep("Variance (\u03b2)", 25),
                       rep("Variance (M)", 25),
                       rep("Variance (\u03b2, Cor)",25),
                       rep("Variance (M, Cor)", 25),
                       rep("PCA",25),
                       rep("KS-like",25),
                       rep("None",25))
)
performanceDF$FeatureSelection <- factor(performanceDF$FeatureSelection,
                                         levels = unique(performanceDF$FeatureSelection))

# Combine into single data frame
performanceDF_test <- data.frame(
  RMSE = unlist(Performances_test),
  FeatureSelection = c(rep("S-score", 1),
                       rep("Variance (\u03b2)", 1),
                       rep("Variance (M)", 1),
                       rep("Variance (\u03b2, Cor)",1),
                       rep("Variance (M, Cor)", 1),
                       rep("PCA",1),
                       rep("KS-like",1),
                       rep("None", 1))
)
performanceDF_test$FeatureSelection <- factor(performanceDF_test$FeatureSelection,
                                         levels = unique(performanceDF_test$FeatureSelection))


  
# Load data
load(paste0("CV_CAIDE1/CV_CAIDE1_Cor_EN.RData"))
load("X/X_Cor_CAIDE1/X_test_Cor.RData")

# Performance in training data
optPar <- which.min(rowMeans(perf))
optPerf <- NULL
for (i in 1:length(trainResults)){
  optPerf <- c(optPerf,trainResults[[i]]$RMSE[optPar])
}
performanceDF <- rbind.data.frame(performanceDF, data.frame(RMSE = optPerf,
                                                            FeatureSelection = rep("Correlation", 25)))

# Performance in test data
testData <- log2(X_test_Cor/(1-X_test_Cor))
pred_test <- predict(finalModel, t(testData))
perf_test <- RMSE(pred = pred_test,
                  obs = Y_test$CAIDE)

performanceDF_test <- rbind.data.frame(performanceDF_test, data.frame(RMSE = perf_test,
                                                            FeatureSelection = rep("Correlation", 1)))

  
performanceDF$FeatureSelection <- factor(performanceDF$FeatureSelection,
                                         levels = c(rep("S-score", 1),
                                         rep("Variance (\u03b2)", 1),
                                         rep("Variance (M)", 1),
                                         rep("Variance (\u03b2, Cor)",1),
                                         rep("Variance (M, Cor)", 1),
                                         rep("PCA",1),
                                         rep("KS-like",1),
                                         rep("Correlation",1),
                                         rep("None", 1)))


# Make plot
p <- ggplot(performanceDF) +
  geom_boxplot(aes(x = FeatureSelection, y = RMSE, fill = FeatureSelection), alpha = 0.3) +
  geom_point(aes(x = FeatureSelection, y = RMSE, color = FeatureSelection), 
             position=position_jitterdodge(jitter.width = 1), size = 2) +
  geom_point(data = performanceDF_test, aes(x = FeatureSelection, y = RMSE), 
             color = "black", size = 5, shape = 18, alpha = 0.7) +
  xlab("Feature Selection Method") +
  ggtitle(Score) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2")[c(1:5,8,6,7)], "red")) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2")[c(1:5,8,6,7)], "red")) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = paste0(Score, "_RMSE_boxplot.png"), width = 10, height = 6)



################################################################################

# Compare performances of feature selection approaches

################################################################################

# Score and feature selection method
Score = "CAIDE1"
methods = c("S", "var", "varM", "varCor", "varMCor", "KS", "PC", "Non")


# Retrieve performance in CV for the different feature selection methods
Performances <- list()
for (i in 1:length(methods)){
  
  # load output
  load(paste0("CV_CAIDE1/CV_", Score, "_", methods[i],".RData"))
  
  # Put performance into list
  Performances[[i]] <- perf

}
names(Performances) <- methods



# Combine into single data frame
performanceDF <- data.frame(
  RMSE = unlist(Performances),
  FeatureSelection = c(rep("S-score", 25),
                       rep("Variance (\u03b2)", 25),
                       rep("Variance (M)", 25),
                       rep("Variance (\u03b2, Cor)",25),
                       rep("Variance (M, Cor)", 25),
                       rep("Kennard-Stone-like", 25),
                       rep("PCA",25),
                       rep("None", 25))
)
performanceDF$FeatureSelection <- factor(performanceDF$FeatureSelection,
                                         levels = unique(performanceDF$FeatureSelection))

load("trainResults_cor.RData")
perf <- matrix(NA, nrow = 1000, ncol = 25)
for (i in 1:length(trainResults)){
  perf[,i] <- trainResults[[i]]$RMSE
}
optPar <- which.min(rowMeans(perf))
optPerf <- NULL
for (i in 1:length(trainResults)){
  optPerf <- c(optPerf,trainResults[[i]]$RMSE[optPar])
}
performanceDF <- rbind.data.frame(performanceDF, data.frame(RMSE = optPerf,
                                                            FeatureSelection = rep("Correlation", 25)))


# Make plot
p <- ggplot(performanceDF) +
  geom_boxplot(aes(x = FeatureSelection, y = RMSE, fill = FeatureSelection), alpha = 0.3) +
  geom_point(aes(x = FeatureSelection, y = RMSE, color = FeatureSelection), 
             position=position_jitterdodge(jitter.width = 1), size = 2) +
  xlab("Feature Selection Method") +
  ggtitle(Score) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"), "red")) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"), "red")) +
  #scale_color_brewer(palette = "Dark2") +
  #scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = paste0(Score, "_RMSE_boxplot.png"), width = 10, height = 6)


################################################################################

# Evaluate on test data

################################################################################

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"

# Load data
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/", f))
}

# Load finalModel
load(paste0("CV_CAIDE1/CV_CAIDE1_", FeatureSelection, ".RData"))

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

testData <- log2(X_test_Cor/(1-X_test_Cor))
test <- predict(finalModel, t(testData))
plot(test, Y_test$CAIDE)
