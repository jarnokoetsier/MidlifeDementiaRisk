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
methods = c("cor", "Cor_spls")
methodNames <- c("ElasticNet", "sPLS")
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
                     Method = rep(methodNames[j], 25))
  
  performanceDF <- rbind.data.frame(performanceDF,temp)
  
  # Performance in test data
  testData <- log2(X_test_Cor/(1-X_test_Cor))
  pred_test <- predict(finalModel, t(testData))
  perf_test <- RMSE(pred = pred_test,
                    obs = Y_test$CAIDE)
  
  temp <- data.frame(RMSE = perf_test,
                     Method = methodNames[j])
  
  performanceDF_test <- rbind.data.frame(performanceDF_test,temp)
  
}




# Make plot
p <- ggplot() +
  geom_boxplot(data = performanceDF, aes(x = Method, y = RMSE, fill = Method), alpha = 0.3) +
  geom_point(data = performanceDF, aes(x = Method, y = RMSE, color = Method), 
             position=position_jitterdodge(jitter.width = 0.3), size = 2) +
  geom_point(data = performanceDF_test, aes(x = Method, y = RMSE), color = "red", size = 5, shape = 18) +
  xlab("") +
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



