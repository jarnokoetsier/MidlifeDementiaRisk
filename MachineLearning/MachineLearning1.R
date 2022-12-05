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
