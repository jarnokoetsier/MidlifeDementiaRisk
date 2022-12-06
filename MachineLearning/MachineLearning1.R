# Clear workspace and console
rm(list = ls())
cat("\014") 

library(glmnet)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
source("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA/MachineLearning/FUN_MachineLearning.R")

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Load Pheno data
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

# Load varM data
files <- list.files('X_varM')
for (f in files){
  load(paste0("X_varM/",f))
}


#*****************************************************************************#
# Model training
#*****************************************************************************#

#=============================================================================#
# FILL IN
#=============================================================================#

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "varM"

# Prepare data (M-values)
X_train = log2(X_CAIDE1_varM/(1-X_CAIDE1_varM))
Y_train = Y_CAIDE1$CAIDE

# Test if samples are in correct order
all(colnames(X_train) == Y_CAIDE1$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5

#=============================================================================#

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           number = nfold, 
                           repeats = nrep, 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = regressionSummary)

# Set grid for lambda (2.5 for CAIDE)
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Use MSE as performance metric
performance_metric = "RMSE"
MLmethod = "glmnet"

# Register cores for parallel computing
#detectCores()
#nCores <- 3
#cl <- makeCluster(nCores)
#registerDoParallel(cl)

# Actual training
set.seed(123)
fit <- train(x = t(X_train),
             y = Y_train,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Stop clusters
#stopCluster(cl)

# Get results
trainResults <- fit$results

# Get optimal lambda and alpha
optAlpha <- fit$bestTune$alpha
optLambda <- fit$bestTune$lambda

# Get coefficient value during the repeated CV
coefs <- matrix(NA, ncol(t(X_train))+1, nfold*nrep)
perf <- rep(NA, nfold*nrep)
folds <- fit$control$index
pred_CV <- NULL
obs_CV <- NULL
fold_CV <- NULL
count = 0
for (r in 1:nrep){
  for (f in 1:nfold){
    count = count + 1
    en_model_cv <- glmnet(x = t(X_train)[folds[[count]],], 
                          y = Y_train[folds[[count]]], 
                          family = "gaussian",
                          alpha = optAlpha, 
                          lambda = optLambda,
                          standardize = TRUE)
    
    # Get coefficients
    coefs[,count] <- as.matrix(coef(en_model_cv))
    
    # Get prediction
    pred <- predict(en_model_cv, t(X_train)[-folds[[count]],])[,1]
    
    # Get performance
    perf[count] <- RMSE(pred = pred, obs = Y_train[-folds[[count]]])
    
    pred_CV <- c(pred_CV,pred)
    obs_CV <- c(obs_CV, Y_train[-folds[[count]]])
    fold_CV <- c(fold_CV, rep(count,length(pred)))
  }
}

ObsPred_CV <- data.frame(Predicted = pred_CV,
                         Observed = obs_CV,
                         Fold = fold_CV)

# Get final model
finalModel <- glmnet(x = t(X_train), 
                     y = Y_train, 
                     family = "gaussian",
                     alpha = optAlpha, 
                     lambda = optLambda,
                     standardize = TRUE)

# Save output
save(trainResults, optLambda, optAlpha, perf, ObsPred_CV, coefs, finalModel,
     file = paste0("CV_", Score, "_", FeatureSelection,".RData"))


#*****************************************************************************#
# Heatmap
#*****************************************************************************#
# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Non"

# load output
load(paste0("CV_", Score, "_", FeatureSelection,".RData"))

# Make heatmap
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
  ggtitle(Score, subtitle = "No Feature Selection") +
  #ggtitle(Score, subtitle = "Variance (M)-based Feature Selection") +
  #ggtitle(Score, subtitle = "Variance (\u03b2)-based Feature Selection") +
  #ggtitle(Score, subtitle = "S-score-based Feature Selection") +
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

ggsave(RMSE_heatmap, file = paste0(Score, "_", FeatureSelection, "_RMSE_heatmap.png"), width = 8, height = 6)


#*****************************************************************************#
# Observed vs predicted
#*****************************************************************************#

ObsVsPred <- ggplot(ObsPred_CV) +
  geom_bin2d(aes(x = Observed, y = Predicted), bins = 100) +
  #ggtitle(Score, subtitle = "S-score-based Feature Selection") +
  #ggtitle(Score, subtitle = "Variance (\u03b2)-based Feature Selection") +
  #ggtitle(Score, subtitle = "Variance (M)-based Feature Selection") +
  ggtitle(Score, subtitle = "No Feature Selection") +
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

ggsave(ObsVsPred, file = paste0(Score, "_", FeatureSelection, "_ObsVsPred.png"), width = 8, height = 6)

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
  ylab("Predicted CAIDE1 Score") +
  ggtitle(Score, subtitle = "No Feature Selection") +
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


R2(ObsPred_CV$Predicted, ObsPred_CV$Observed)

#*****************************************************************************#
# Regression coefficients
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
load("X_coefs.RData")

# Perform PCA
pcaList <-  prcomp(t(X_coefs),        
                   retx = TRUE,
                   center =TRUE,
                   scale = TRUE,
                   rank. = 10)

# Save pcaList object
#save(pcaList, file = "pcaList_S.RData")
#load("pcaList_S.RData")

# Get PCA scores
PCAscores <- as.data.frame(pcaList$x)
PCAscores$ID <- rownames(PCAscores)

# Combine with sample info
PCAscores <- inner_join(PCAscores, dat, by = c("ID" = "Basename"))

# Combine with cell type composition
PCAscores <- inner_join(PCAscores,cellType, by = c("ID" = "ID"))

# Combine with cAIDE
load("E:/Thesis/EXTEND/Phenotypes/CAIDE.Rdata")
PCAscores <- inner_join(PCAscores,CAIDE, by = c("ID" = "Basename"))

# Get explained variance
explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)



p_12 <- ggplot(data = PCAscores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.9, size = 2, aes(color = Age.x)) +
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



# Get Principal Components
PCs <- PCAscores[,1:5]
colnames(PCs) <- paste0(colnames(PCs), " (", explVar[1:5], "%)")

# Get sample information
meta_cell <- PCAscores[,colnames(cellType)[-c(7,8)]]
meta_caide <- PCAscores[,176:183]

colnames(meta_cell) <- c("CD8 T-cells", "CD4 T-cells", "NK cells", "B-cells", "Monocytes",
                    "Neutrophils")

colnames(meta_caide) <- c("Age", "Sex", "Education", "Blood Pressure", "BMI",
                         "Serum HDL", "Physical Activity", "CAIDE1")

meta <- cbind.data.frame(meta_caide, meta_cell)

# Get correlations
corDF <- expand.grid(colnames(PCs), colnames(meta))
colnames(corDF) <- c("PCs", "Meta")
corDF$Correlation <- rep(NA, nrow(corDF))
for (i in 1:nrow(corDF)){
  corDF$Correlation[i] <- cor(PCs[,corDF$PCs[i]], meta[,corDF$Meta[i]], method = "spearman")
}

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

ggsave(p,file = "correlationPlot_coefsModel.png", width = 8, height = 6)


#*****************************************************************************#
# Multiple regression
#*****************************************************************************#
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




# Multiple regression  (CAIDE factors)
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
#*****************************************************************************#
# Performance in CV
#*****************************************************************************#

#=============================================================================#
# CAIDE1
#=============================================================================#

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "S"

# load output
load(paste0("CV_", Score, "_", FeatureSelection,".RData"))

# Get performance
CAIDE1_S_non_perf <- perf

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "var"

# load output
load(paste0("CV_", Score, "_", FeatureSelection,".RData"))

# Get performance
CAIDE1_var_non_perf <- perf

# Combine into single data frame
performanceDF_CAIDE1 <- data.frame(
  RMSE = c(CAIDE1_S_non_perf, CAIDE1_var_non_perf),
  FeatureSelection = c(rep("S-score", 25),
                       rep("Variance", 25)),
  Correction = c(rep("No Correction", 25),
                 rep("No Correction", 25))
)

ggplot(performanceDF_CAIDE1) +
  geom_boxplot(aes(x = FeatureSelection, y = RMSE, fill = Correction), alpha = 0.3) +
  geom_point(aes(x = FeatureSelection, y = RMSE, color = Correction), 
             position=position_jitterdodge(jitter.width = 0.2), size = 2) +
  xlab("Feature Selection Method") +
  ggtitle("CAIDE1") +
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
