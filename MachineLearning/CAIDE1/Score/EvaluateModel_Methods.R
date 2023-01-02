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

