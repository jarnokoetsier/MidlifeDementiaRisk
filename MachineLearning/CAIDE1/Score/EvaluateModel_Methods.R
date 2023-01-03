# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(kernlab)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(patchwork)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load cell type composition and phenotype data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")



################################################################################

# Compare performances of methods

################################################################################

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"

# Load data
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/", f))
}
# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}


# Load finalModel
performanceDF <- NULL
performanceDF_test <- NULL
methods = c("cor", "Cor_spls", "Cor_svm")
methodNames <- c("ElasticNet", "sPLS", "SVM")
selectionNames <- c("Correlation", "Correlation", "Correlation")
for (j in 1:length(methods)){
  
  # Load data
  load(paste0("CV_CAIDE1/CV_CAIDE1_", methods[j], ".RData"))
  
  # Performance in training data
  optPar <- which.min(rowMeans(perf))
  optPerf <- NULL
  for (i in 1:length(trainResults)){
    optPerf <- c(optPerf,trainResults[[i]]$RMSE[optPar])
  }
  temp <- data.frame(RMSE = optPerf,
                     Method = rep(methodNames[j], 25),
                     Selection = rep(selectionNames[j],25))
  
  performanceDF <- rbind.data.frame(performanceDF,temp)
  
  # Performance in test data
  testData <- log2(X_test_Cor/(1-X_test_Cor))
  pred_test <- predict(finalModel, t(testData))
  perf_test <- RMSE(pred = pred_test,
                    obs = Y_test$CAIDE)
  
  temp <- data.frame(RMSE = perf_test,
                     Method = methodNames[j],
                     Selection = selectionNames[j])
  
  performanceDF_test <- rbind.data.frame(performanceDF_test,temp)
  
}




# Make plot
p <- ggplot() +
  geom_boxplot(data = performanceDF, aes(x = Method, y = RMSE, fill = Method), alpha = 0.3) +
  geom_point(data = performanceDF, aes(x = Method, y = RMSE, color = Method), 
             position=position_jitterdodge(jitter.width = 0.3), size = 2) +
  geom_point(data = performanceDF_test, aes(x = Method, y = RMSE), 
             color = "black", size = 5, shape = 18, alpha = 0.7) +
  xlab("") +
  ggtitle(Score) +
  #scale_fill_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"), "red")) +
  #scale_color_manual(values = c(RColorBrewer::brewer.pal(8, "Dark2"), "red")) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = paste0(Score, "_RMSE_boxplot_methods.png"), width = 10, height = 6)


################################################################################

# Compare features

################################################################################

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"

# Load data
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/", f))
}
# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Score and feature selection method

methods <- c("cor", "Cor_spls")
methodNames <- c("ElasticNet", "sPLS")
selectionNames <- c("Correlation", "Correlation")
coefficients <- matrix(NA,nrow = 10000,ncol=2)
rownames(coefficients) <- rownames(X_CAIDE1_Cor)
for (j in 1:length(methods)){
  
  # Load data
  load(paste0("CV_CAIDE1/CV_CAIDE1_", methods[j], ".RData"))
  
  # Get coefficients
  coefficients[,j] <- coef(finalModel)[rownames(coefficients),1]
  
}
colnames(coefficients) <- methodNames
coefficients <- as.data.frame(coefficients)

# variable importance svm



p <- ggplot(coefficients) +
  geom_point(aes(x = ElasticNet, y = sPLS), alpha = 0.5) +
  scale_color_viridis_c() +
  xlab("Regression coefficient (ElasticNet)") +
  ylab("Regression coefficient (sPLS)") +
  theme_minimal()

ggsave(p, file = "varImp_CorModels.png",width = 8, height = 8)



# TEST #
dataMatrix <- cbind(finalModel$projection, test[rownames(finalModel$projection),1])
model <- lm(V10 ~0 + .,data = as.data.frame(dataMatrix))

load("CV_CAIDE1/CV_CAIDE1_Cor_svm.RData")
test <- varImp(finalModel)

################################################################################

# sPLS components

################################################################################

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"

load(paste0("CV_CAIDE1/CV_CAIDE1_", "Cor_spls", ".RData"))

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Calculate scores
scores <- finalModel$x[,rownames(finalModel$projection)] %*% finalModel$projection
colnames(scores) <- paste0("Component ", 1:9)
all(rownames(scores) == Y_CAIDE1$Basename)

# Get meta data
meta <- Y_CAIDE1[,14:21]
colnames(meta) <- c("Age", "Sex", "Education", "Blood Pressure", 
                    "BMI", "Cholesterol", "Physical", "CAIDE1")

# Get correlations
corDF <- expand.grid(colnames(scores), colnames(meta))
colnames(corDF) <- c("PCs", "Meta")
corDF$Correlation <- rep(NA, nrow(corDF))
for (i in 1:nrow(corDF)){
  corDF$Correlation[i] <- cor(scores[,corDF$PCs[i]], meta[,corDF$Meta[i]], method = "spearman")
}

# Make plot
p <- ggplot() +
  geom_point(data = corDF, aes(x = PCs, y = Meta, color = Correlation, size = abs(Correlation))) +
  labs(color = "Spearman\nCorrelation", size = "|Spearman\nCorrelation|") +
  guides(size = "none") +
  #ggtitle(subtitle) +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))

# save plot
ggsave(p,file = "correlationPlot_sPLS.png", width = 8, height = 6)


################################################################################

# Cohen's f2 (Elastic Net)

################################################################################


load(paste0("CV_CAIDE1/CV_CAIDE1_", "cor", ".RData"))
load(paste0("X/X_Cor/X_CAIDE1_Cor.RData"))

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Get top 50 features
topFeatures <- names(tail(sort(abs(coef(finalModel)[-1,1])),50))
X_CAIDE1_top <- X_CAIDE1_Cor[topFeatures,]
topCoefs <- data.frame(CpG = topFeatures,
                         coefValue = coef(finalModel)[topFeatures,1])

# Check whether samples are in correct order
all(colnames(X_CAIDE1_top) == Y_CAIDE1$Basename)

# Make formula
formula <- paste0("cbind(",paste(topFeatures, collapse = ", "),") ~ ", 
                  paste0("0 + ", paste(colnames(Y_CAIDE1[,14:20]), collapse = " + ")))

# Mean center the data
X_coefs_scaled <- (X_CAIDE1_top - rowMeans(X_CAIDE1_top))

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
  formula <- paste0("cbind(",paste(topFeatures, collapse = ", "),") ~ ", 
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
factorNames <- c("Age", "Sex", "Edu", "BP", "BMI", "Chol", "Physical")
plotDF <- data.frame(cohenF = c(unlist(cohenF)/globalEffect, globalEffect),
                     Effect = rep(c(factorNames, "Global"), each = nrow(X_coefs_scaled)),
                     CpG = rep(topFeatures,length(factorNames) +1))

# Reorder
plotDF$Effect <- factor(plotDF$Effect, levels = c(factorNames, "Global"))

# Combine with coefficient values in final model
plotDF <- inner_join(plotDF, topCoefs, by = c("CpG" = "CpG"))
plotDF$CpG <- factor(plotDF$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))


# Make plot
main <- ggplot(plotDF[plotDF$Effect != "Global",]) +
  geom_bar(aes(y = cohenF, x = CpG, fill = Effect), stat="identity", alpha = 1) +
  facet_grid(rows = vars(Effect)) +
  xlab(NULL) +
  ylab(expression("Cohen's " ~ f^2 ~ " (local / global)")) +
  ggtitle("") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill= "grey", linewidth = 1, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 13))

global <- ggplot(plotDF[plotDF$Effect == "Global",]) +
  geom_bar(aes(y = cohenF, x = CpG, fill = Effect), stat="identity") +
  xlab(NULL) +
  ylab(expression("Cohen's " ~ f^2)) +
  ggtitle("") +
  facet_grid(rows = vars(Effect)) +
  #scale_fill_manual(values = "black") +
  scale_fill_manual(values = "#666666") +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.15,0.2),labels = c("0.0", "", "0.1", "", "0.2")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        strip.background = element_rect(fill= "#ff726f",
                                        linewidth = 1, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 13))



topCoefs$CpG <- factor(topCoefs$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))
top <- ggplot(topCoefs) +
  #geom_point(aes(x = fct_reorder(cpg, avgValue), y = avgValue), color = "blue") +
  geom_bar(aes(x = CpG, y = coefValue, fill = coefValue), stat = "identity", color = "black") +
  ylab("Coefficients\nFinal Model") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.5,0.5), oob = scales::squish) +
  #scale_color_viridis_c(limits = c(-0.5, 0.5), oob = scales::squish) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p <- top / main / global +
  plot_layout(heights = c(1,7,1))

p <- main / global +
  plot_layout(heights = c(7,1))

# Save plot
ggsave(p,file = "cohenF_ElasticNetModel_Cor_CAIDE1.png", width = 8, height = 8)


################################################################################

# Cohen's f2 (sPLS)

################################################################################


load(paste0("CV_CAIDE1/CV_CAIDE1_", "Cor_spls", ".RData"))
load(paste0("X/X_Cor/X_CAIDE1_Cor.RData"))

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

# Get top 50 features
topFeatures <- names(tail(sort(abs(coef(finalModel)[-1,1])),50))
X_CAIDE1_top <- X_CAIDE1_Cor[topFeatures,]
topCoefs <- data.frame(CpG = topFeatures,
                       coefValue = coef(finalModel)[topFeatures,1])

# Check whether samples are in correct order
all(colnames(X_CAIDE1_top) == Y_CAIDE1$Basename)

# Make formula
formula <- paste0("cbind(",paste(topFeatures, collapse = ", "),") ~ ", 
                  paste0("0 + ", paste(colnames(Y_CAIDE1[,14:20]), collapse = " + ")))

# Mean center the data
X_coefs_scaled <- (X_CAIDE1_top - rowMeans(X_CAIDE1_top))

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
  formula <- paste0("cbind(",paste(topFeatures, collapse = ", "),") ~ ", 
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
factorNames <- c("Age", "Sex", "Edu", "BP", "BMI", "Chol", "Physical")
plotDF <- data.frame(cohenF = c(unlist(cohenF)/globalEffect, globalEffect),
                     Effect = rep(c(factorNames, "Global"), each = nrow(X_coefs_scaled)),
                     CpG = rep(topFeatures,length(factorNames) +1))

# Reorder
plotDF$Effect <- factor(plotDF$Effect, levels = c(factorNames, "Global"))

# Combine with coefficient values in final model
plotDF <- inner_join(plotDF, topCoefs, by = c("CpG" = "CpG"))
plotDF$CpG <- factor(plotDF$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))


# Make plot
main <- ggplot(plotDF[plotDF$Effect != "Global",]) +
  geom_bar(aes(y = cohenF, x = CpG, fill = Effect), stat="identity", alpha = 1) +
  facet_grid(rows = vars(Effect)) +
  xlab(NULL) +
  ylab(expression("Cohen's " ~ f^2 ~ " (local / global)")) +
  ggtitle("") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill= "grey", linewidth = 1, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 13))

global <- ggplot(plotDF[plotDF$Effect == "Global",]) +
  geom_bar(aes(y = cohenF, x = CpG, fill = Effect), stat="identity") +
  xlab(NULL) +
  ylab(expression("Cohen's " ~ f^2)) +
  ggtitle("") +
  facet_grid(rows = vars(Effect)) +
  #scale_fill_manual(values = "black") +
  scale_fill_manual(values = "#666666") +
  #scale_y_continuous(breaks = c(0,0.05,0.1,0.15,0.2),labels = c("0.0", "", "0.1", "", "0.2")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        strip.background = element_rect(fill= "#ff726f",
                                        linewidth = 1, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 13))



topCoefs$CpG <- factor(topCoefs$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))
top <- ggplot(topCoefs) +
  #geom_point(aes(x = fct_reorder(cpg, avgValue), y = avgValue), color = "blue") +
  geom_bar(aes(x = CpG, y = coefValue, fill = coefValue), stat = "identity", color = "black") +
  ylab("Coefficients\nFinal Model") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.03,0.03), oob = scales::squish) +
  #scale_color_viridis_c(limits = c(-0.5, 0.5), oob = scales::squish) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

p <- top / main / global +
  plot_layout(heights = c(1,7,1))

p <- main / global +
  plot_layout(heights = c(7,1))

# Save plot
ggsave(p,file = "cohenF_sPLStModel_Cor_CAIDE1.png", width = 8, height = 8)

