# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")


m <- c("KS", "Cor", "CorKS", "CorKS1")
methods <- c("KS-like","Correlation", "Corr (20,000) + KS-like", "Corr (40,000) + KS-like")

plotExplVar = NULL
for (i in 1:length(methods)){
  
  files <- list.files(paste0("X/X_", m[i]))
  for (f in files){
    load(paste0("X/X_", m[i], "/",f))
  }
  
  # Get data
  if (methods[i] == "KS-like"){
    dataMatrix <- X_nonTest_KS
  }
  if (methods[i] == "Correlation"){
    dataMatrix <- X_nonTest_Cor
  }
  if (methods[i] == "Corr (20,000) + KS-like"){
    dataMatrix <- X_nonTest_CorKS
  }
  if (methods[i] == "Corr (40,000) + KS-like"){
    dataMatrix <- X_nonTest_CorKS
  }
  
  # Perform PCA
  pcaList <-  prcomp(t(dataMatrix),        
                     retx = TRUE,
                     center =TRUE,
                     scale = TRUE)
  
  # Get explained variance
  explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)
  cumVar <- sapply(1:length(explVar),function(x){sum(explVar[1:x])})
  
  temp <- data.frame(
    explVar = explVar,
    cumVar = cumVar,
    FeatureSelection = rep(methods[i], length(explVar))
  )
  
  plotExplVar <- rbind.data.frame(plotExplVar, temp)
}
plotExplVar$PC <- factor(rep(paste0("PC", 1:924), length(methods)),
                         levels = paste0("PC", 1:924))

plotExplVar$FeatureSelection <- factor(plotExplVar$FeatureSelection,
                                       levels = methods)

color <- c(RColorBrewer::brewer.pal(n = 8, "Dark2")[7],"#FB5607","#C1121F")
p  <- ggplot(plotExplVar[plotExplVar$FeatureSelection != "KS-like",]) +
  geom_step(aes(x = PC, y = cumVar, group = FeatureSelection, color = FeatureSelection),
            linewidth = 1.5) +
  geom_hline(yintercept = 100,  linetype = "dashed", linewidth = 1.5) +
  xlab("Principal Components (PC1 - 924)") +
  ylab("Cumulative Explained Variance (%)") +
  #scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = color) +
  ylim(c(0,100)) +
  theme_classic() +
  theme(legend.title  = element_blank(),
        legend.position = "bottom",
        axis.text.x = element_blank())#element_text(angle = 45, vjust = 0.5))

ggsave(p, file = "FeatureSelectionExplVar_CorKS1.png", width = 6, height = 5)


###############################################################################

# 3. Venn diagram

###############################################################################

library(ggvenn)

m1 <- "CorKS"
files <- list.files(paste0("X/X_", m1))
for (f in files){
  load(paste0("X/X_", m1, "/",f))
}
plotList <- list(rownames(X_nonTest_Cor), 
                 rownames(X_nonTest_CorKS), 
                 rownames(dataMatrix))

names(plotList) <- c("Correlation","Corr (20,000) + KS-like", "Corr (40,000) + KS-like")


p <- ggvenn(plotList, show_percentage = FALSE, fill_color = color)
ggsave(p, file = "FeatureSelectionVenn_CorKS.png", width = 6, height = 6)


###############################################################################

# Performance CAIDE1 prediction

###############################################################################
library(glmnet)
library(spls)
library(kernlab)
library(caret)

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}


# Load finalModel
performanceDF <- NULL
performanceDF_test <- NULL
methods = c("cor", "Cor_spls", "Cor_svm", 
            "CorKS", "CorKS_spls", "CorKS_svm", 
            "CorKS1", "CorKS1_spls", "CorKS1_svm")
methodNames <- rep(c("ElasticNet", "sPLS", "SVM"),3)
selectionNames <- c(rep("Correlation",3), 
                    rep("Corr (20,000) + KS-like",3), 
                    rep("Corr (40,000) + KS-like",3))
for (j in 1:length(methods)){
  
  # Load data
  file <- paste0("CV_CAIDE1/CV_CAIDE1_", methods[j], ".RData")
  
  if (file.exists(file)){
    load(file)
    
    if (selectionNames[j] == "Correlation"){
      
      FeatureSelection <- "Cor"
      files <- list.files(paste0("X/X_", FeatureSelection))
      for (f in files){
        load(paste0("X/X_", FeatureSelection, "/",f))
      }
      
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
    if (selectionNames[j] == "Corr (20,000) + KS-like"){
     
      FeatureSelection <- "CorKS"
      files <- list.files(paste0("X/X_", FeatureSelection))
      for (f in files){
        load(paste0("X/X_", FeatureSelection, "/",f))
      }
      
      
      # Performance in training data
      temp <- data.frame(RMSE = perf,
                         Method = rep(methodNames[j], 25),
                         Selection = rep(selectionNames[j],25))
      
      performanceDF <- rbind.data.frame(performanceDF,temp)
      
      # Performance in test data
      testData <- log2(X_test_CorKS/(1-X_test_CorKS))
      pred_test <- predict(finalModel, t(testData))
      perf_test <- RMSE(pred = pred_test,
                         obs = Y_test$CAIDE)
      
      temp <- data.frame(RMSE = perf_test,
                         Method = methodNames[j],
                         Selection = selectionNames[j])
      
      performanceDF_test <- rbind.data.frame(performanceDF_test,temp)
    }
    if (selectionNames[j] == "Corr (40,000) + KS-like"){
      
      FeatureSelection <- "CorKS1"
      files <- list.files(paste0("X/X_", FeatureSelection))
      for (f in files){
        load(paste0("X/X_", FeatureSelection, "/",f))
      }
      
      # Performance in training data
      temp <- data.frame(RMSE = perf,
                         Method = rep(methodNames[j], 25),
                         Selection = rep(selectionNames[j],25))
      
      performanceDF <- rbind.data.frame(performanceDF,temp)
      
      # Performance in test data
      testData <- log2(X_test_CorKS/(1-X_test_CorKS))
      pred_test <- predict(finalModel, t(testData))
      perf_test <- RMSE(pred = pred_test,
                        obs = Y_test$CAIDE)
      
      temp <- data.frame(RMSE = perf_test,
                         Method = methodNames[j],
                         Selection = selectionNames[j])
      
      performanceDF_test <- rbind.data.frame(performanceDF_test,temp)
    }
    
  }
}




# Make plot
performanceDF$Selection <- factor(performanceDF$Selection, 
                                  levels = c("Correlation", "Corr (20,000) + KS-like", "Corr (40,000) + KS-like"))
performanceDF_test$Selection <- factor(performanceDF_test$Selection, 
                                  levels = c("Correlation", "Corr (20,000) + KS-like", "Corr (40,000) + KS-like"))

color <- c(RColorBrewer::brewer.pal(n = 8, "Dark2")[7],"#FB5607","#C1121F")
p <- ggplot() +
  geom_boxplot(data = performanceDF, aes(x = Method, y = RMSE, fill = Selection), alpha = 0.3) +
  geom_point(data = performanceDF, aes(x = Method, y = RMSE, color = Selection), 
             position=position_jitterdodge(jitter.width = 0.1), size = 2) +
  geom_point(data = performanceDF_test, aes(x = Method, y = RMSE, shape = Selection), 
             position=position_jitterdodge(jitter.width = 0), color = "black", size = 5, alpha = 0.7) +
  xlab("") +
  ggtitle("CAIDE1") +
  scale_fill_manual(values = color) +
  scale_color_manual(values = color) +
  #scale_color_brewer(palette = "Set1") +
  #scale_fill_brewer(palette = "Set1") +
  scale_shape_manual(values = c(18,18,18), guide = "none")+
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = paste0("CAIDE1_RMSE_boxplot_methods.png"), width = 10, height = 6)

