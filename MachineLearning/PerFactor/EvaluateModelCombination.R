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
library(RColorBrewer)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Meta data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

load("PerFactor/CombineFactors/predictedScore_factors_EXTEND.RData")
load(paste0("PerFactor/CombineFactors/Fit_CombineFactors_", "CAIDE1", "_", "RF", ".RData"))


# make prediction
X_test <- predictedScore_factors[Y_test$Basename,]
pred <- predict(fit, X_test)

plotDF <- data.frame(Predicted = pred,
                     Observed = Y_test$CAIDE)

R2(obs = plotDF$Observed, pred = plotDF$Predicted)

p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(x = Predicted, y = Observed, color = Observed-Predicted), size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("CAIDE1") +
  xlab("Observed Score") +
  ylab("Predicted Score") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 


ggsave(p,file = "ObsvsPred_DataDrivenCombination.png", width = 8, height = 4)