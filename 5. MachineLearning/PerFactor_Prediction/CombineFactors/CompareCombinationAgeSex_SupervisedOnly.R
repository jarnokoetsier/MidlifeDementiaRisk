# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(caret)
library(glmnet)
library(spls)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Make vectors to save performances
perf_CAIDE1 <- rep(NA,3)
perf_CAIDE2 <- rep(NA,3)
perf_LIBRA <- rep(NA,3)

# Calculate performances in test set
methods <- c("EN", "sPLS", "RF")
load("Y/Y_test.RData")
load("PerFactor/CombineFactors/predictedScore_factors_EXTEND.RData")
for (i in 1:length(methods)){
  
  # CAIDE1
  
  # Load model
  load(paste0("PerFactor/CombineFactors/Fit_CombineFactors_", "CAIDE1", "_", methods[i], ".RData"))
  
  # Test set
  X_test <- predictedScore_factors[Y_test$Basename,]
  
  # Set age and sex to mean value
  X_test1 <- X_test
  X_test1$Age <- rep(mean(X_test$Age),nrow(X_test))
  X_test1$SexMale <- rep(mean(X_test$SexMale),nrow(X_test))
  pred <- predict(fit, X_test1)
  perf_CAIDE1[i] <- R2(pred = pred, obs = Y_test$CAIDE)
  
  # CAIDE2
  
  # Load model
  load(paste0("PerFactor/CombineFactors/Fit_CombineFactors_", "CAIDE2", "_", methods[i], ".RData"))
  
  # make test set
  X_test <- predictedScore_factors[Y_test$Basename,]
  
  # Set age and sex to mean value
  X_test1 <- X_test
  X_test1$Age <- rep(mean(X_test$Age),nrow(X_test))
  X_test1$SexMale <- rep(mean(X_test$SexMale),nrow(X_test))
  pred <- predict(fit, X_test1)
  perf_CAIDE2[i] <- R2(pred = pred, obs = Y_test$CAIDE2)
  
  # LIBRA
  
  # Load model
  load(paste0("PerFactor/CombineFactors/Fit_CombineFactors_", "LIBRA", "_", methods[i], ".RData"))
  
  # make test set
  X_test <- predictedScore_factors[Y_test$Basename,]
  
  # Set age and sex to mean value
  X_test1 <- X_test
  X_test1$Age <- rep(mean(X_test$Age),nrow(X_test))
  X_test1$SexMale <- rep(mean(X_test$SexMale),nrow(X_test))
  pred <- predict(fit, X_test1)
  perf_LIBRA[i] <- R2(pred = pred, obs = Y_test$LIBRA)
  
}

# Prepare data for plotting
plotDF <- data.frame(R2 = c(perf_CAIDE1, perf_CAIDE2, perf_LIBRA),
                     Method = rep(c("ElasticNet", "sPLS", "Random Forest")),
                     Score = c(rep("CAIDE1",3), rep("CAIDE2",3), rep("LIBRA",3)))

plotDF$Method <- factor(plotDF$Method, 
                        levels = rev(c(
                                       "ElasticNet",
                                       "sPLS", 
                                       "Random Forest"
                        )))

# Set colors
colors <- c(RColorBrewer::brewer.pal(n = 8, name = "Reds")[4:8],
            RColorBrewer::brewer.pal(n = 8, name = "Oranges")[4:8],
            RColorBrewer::brewer.pal(n = 8, name = "PuRd")[4:8])

# Make plot
p <- ggplot(plotDF[(plotDF$Method == "Random Forest") |
                     (plotDF$Method == "sPLS") |
                     (plotDF$Method == "ElasticNet"), ]) +
  geom_bar(aes(x = Method, y = R2, fill = Score, alpha = Method), stat = "identity", color = "black") +
  facet_grid(rows = vars(Score), scale = "free", space = "free") +
  coord_flip() +
  theme_bw() +
  xlab(NULL) +
  ylab(expression(R^2)) +
  ylim(c(0,0.5)) +
  scale_alpha_manual(values = rev(c(0.6,0.8,1))) +
  scale_fill_manual(values = c("#EF3B2C","#FE9929", "#807DBA")) +
  theme(legend.position = "none")

# Save plot
ggsave(p, file = "CombinationMethods_supervised_noAgeSex.png", width = 8, height = 6)