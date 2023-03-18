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
  geom_boxplot(aes(x = Method, y = R2, fill = Score, alpha = Feature), 
               outlier.shape = NA) +
  geom_point(aes(x = Method, y = R2, color = Score, shape = Feature), 
             position=position_jitterdodge(jitter.width = 0.1), 
             size = 2, alpha = 0.8) +
  geom_point(data = plotDF_test_all, aes(x = Method, y = R2, shape = Feature), 
             fill = "black", position=position_jitterdodge(jitter.width = 0), 
             size = 4, alpha = 0.7) +
  facet_grid(cols = vars(Score)) +
  ylab(expression(R^2)) +
  xlab(NULL) +
  guides(color = "none", fill = "none") +
  scale_shape_manual(values = c(18,20)) +
  scale_alpha_manual(values = c(0.5, 1)) +
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


ggsave(p,file = "Evaluate_CorPGS1.png", width = 8, height = 5)




plotDF_test_all$Group <- paste0(plotDF_test_all$Score, "_", plotDF_test_all$Method)

testDF <- plotDF_test_all %>%
  group_by(Group) %>%
  summarise(Color = -1*sign(diff(R2)),
            Method = Method,
            Score = Score,
            R2 = R2,
            Feature = Feature)

testDF$Color <- as.character(testDF$Color)
testDF$Color1 <- paste0(testDF$Color,"_",testDF$Feature)

p <- ggplot(testDF) +
  geom_line(aes(x = R2, y = Method, group = Group, color = Color), linewidth = 1.5) +
  geom_point(aes(x = R2,y = Method, shape = Feature, size = Feature, color = Color1)) +
  facet_grid(rows = vars(Score)) +
  guides(color = "none") +
  xlab(expression(R^2)) +
  ylab(NULL)+
  scale_size_manual(values = c(3,3)) +
  scale_shape_manual(values = c(15,16)) +
  scale_color_manual(values = c("red","black", "red", "blue","black", "blue")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave(p,file = "Evaluate_CorPGS_test.png", width = 6, height = 6)


  
  