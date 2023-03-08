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

# Combine with PGS
load("E:/Thesis/EXTEND/df_list.RData")
PGS <- df_list$bayesr
samples_test <- Y_test[Y_test$ID %in% rownames(PGS),]
PGS_test <- t(PGS[samples_test$ID,])
colnames(PGS_test) <- samples_test$Basename


###############################################################################

# All scores

###############################################################################

methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")
scores <- c("CAIDE1", "CAIDE2", "LIBRA")
plotDF <- NULL
plotDF_test <- NULL
for (s in 1:length(scores)){
  for (m in 1:length(methods)){
    # Load data
    if (file.exists(paste0("Hybrid/CV_",scores[s],"_CorPGS_", methods[m],".RData"))){
      load(paste0("Hybrid/CV_",scores[s],"_CorPGS_", methods[m],".RData"))
      
      if (scores[s] == "CAIDE1"){
        load("X/X_Cor_CAIDE1/X_test_Cor.RData")
        X_test <- log2(X_test_Cor/(1-X_test_Cor))
        testData <- rbind(X_test[,samples_test$Basename], PGS_test)
        test <- predict(finalModel,t(testData))
        
        perf_test <- R2(pred = test, obs = samples_test$CAIDE)
      }
      if (scores[s] == "CAIDE2"){
        load("X/X_Cor_CAIDE2/X_test_Cor2.RData")
        X_test <- log2(X_test_Cor2/(1-X_test_Cor2))
        testData <- rbind(X_test[,samples_test$Basename], PGS_test)
        test <- predict(finalModel,t(testData))
        
        perf_test <- R2(pred = test, obs = samples_test$CAIDE2)
      }
      if (scores[s] == "LIBRA"){
        load("X/X_Cor_LIBRA/X_test_CorL.RData")
        X_test <- log2(X_test_CorL/(1-X_test_CorL))
        testData <- rbind(X_test[,samples_test$Basename], PGS_test)
        test <- predict(finalModel,t(testData))
        
        perf_test <- R2(pred = test, obs = samples_test$LIBRA)
      }
      temp_test <- data.frame(R2 = perf_test,
                              Score = scores[s],
                              Method = methodNames[m])
      plotDF_test <- rbind.data.frame(plotDF_test, temp_test)
      
      
      # Performance in CV
      optPar <- which.min(rowMeans(perf))
      optPerf <- NULL
      for (i in 1:length(trainResults)){
        optPerf <- c(optPerf,trainResults[[i]]$Rsquared[optPar])
      }
      
      temp <- data.frame(R2 = optPerf,
                         Score = rep(scores[s], length(optPerf)),
                         Method = rep(methodNames[m],length(optPerf)))
      plotDF <- rbind.data.frame(plotDF, temp)
    }
    
  }
  
}


plotDF$Method <- factor(plotDF$Method,
                        levels = methodNames)
plotDF_test$Method <- factor(plotDF_test$Method,
                             levels = methodNames)


plotDF_test$Feature <- rep("CpG + PGS", nrow(plotDF_test))
plotDF$Feature <- rep("CpG + PGS", nrow(plotDF))



methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")
scores <- c("CAIDE1", "CAIDE2", "LIBRA")
plotDF1 <- NULL
plotDF_test1 <- NULL
for (s in 1:length(scores)){
  for (m in 1:length(methods)){
    # Load data
    if (file.exists(paste0("CV_",scores[s],"/CV_",scores[s],"_Cor_", methods[m],".RData"))){
      load(paste0("CV_",scores[s],"/CV_",scores[s],"_Cor_", methods[m],".RData"))
      load(paste0("CVindex_",scores[s],".RData"))
      
      if (scores[s] == "CAIDE1"){
        load("X/X_Cor_CAIDE1/X_test_Cor.RData")
        testData <- log2(X_test_Cor/(1-X_test_Cor))
        pred_test <- predict(finalModel, t(testData))
        perf_test <- R2(pred = pred_test,
                        obs = Y_test$CAIDE)
      }
      if (scores[s] == "CAIDE2"){
        load("X/X_Cor_CAIDE2/X_test_Cor2.RData")
        testData <- log2(X_test_Cor2/(1-X_test_Cor2))
        pred_test <- predict(finalModel, t(testData))
        perf_test <- R2(pred = pred_test,
                        obs = Y_test$CAIDE2)
      }
      if (scores[s] == "LIBRA"){
        load("X/X_Cor_LIBRA/X_test_CorL.RData")
        testData <- log2(X_test_CorL/(1-X_test_CorL))
        pred_test <- predict(finalModel, t(testData))
        perf_test <- R2(pred = pred_test,
                        obs = Y_test$LIBRA)
      }
      temp_test <- data.frame(R2 = perf_test,
                              Score = scores[s],
                              Method = methodNames[m])
      plotDF_test1 <- rbind.data.frame(plotDF_test1, temp_test)
      
      
      # Performance in CV
      optPar <- which.min(rowMeans(perf))
      optPerf <- NULL
      for (i in 1:length(trainResults)){
        optPerf <- c(optPerf,trainResults[[i]]$Rsquared[optPar])
      }
      
      temp <- data.frame(R2 = optPerf,
                         Score = rep(scores[s], length(optPerf)),
                         Method = rep(methodNames[m],length(optPerf)))
      plotDF1 <- rbind.data.frame(plotDF1, temp)
    }
    
  }
  
}


plotDF1$Method <- factor(plotDF1$Method,
                        levels = methodNames)
plotDF_test1$Method <- factor(plotDF_test1$Method,
                             levels = methodNames)

plotDF_test1$Feature <- rep("CpG", nrow(plotDF_test1))
plotDF1$Feature <- rep("CpG", nrow(plotDF1))


plotDF_all <- rbind.data.frame(plotDF, plotDF1)
plotDF_test_all <- rbind.data.frame(plotDF_test, plotDF_test1)

p <- ggplot(plotDF_all) +
  geom_boxplot(aes(x = Feature, y = R2, fill = Method), 
               alpha = 1, outlier.shape = NA) +
  geom_point(aes(x = Feature, y = R2, color = Method), 
             position=position_jitterdodge(jitter.width = 0.1), 
             size = 2, alpha = 0.8) +
  geom_point(data = plotDF_test_all, aes(x = Feature, y = R2, color = Method), 
             fill = "black", position=position_jitterdodge(jitter.width = 0), 
             size = 4, shape = 23, alpha = 0.7) +
  facet_grid(cols = vars(Score)) +
  ylab(expression(R^2)) +
  xlab(NULL) +
  guides(color = "none") +
  scale_fill_manual(values = c("#FB6A4A","#FEC44F", "#BCBDDC")) +
  scale_color_manual(values = c("#A50F15","#CC4C02", "#6A51A3")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

g <- ggplot_gtable(ggplot_build(p))

strips <- which(grepl('strip-', g$layout$name))

pal <- c("#A50F15","#CC4C02", "#6A51A3")


for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- pal[i]
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "white"
}

plot(g)


ggsave(p,file = "RiskScores_test.png", width = 8, height = 5)