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
library(corrr)
library(patchwork)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Meta data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

#==============================================================================#
# CAIDE1
#==============================================================================#
# Get test data
load("X/X_Cor_CAIDE1/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))

accuracy <- rep(NA,6)
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

for (i in 1:length(methods)){
  # Load model
  load(paste0("CV_CAIDE1/CV_CAIDE1_Cor_", methods[i],".RData"))
  
  # make prediction
  pred <- predict(finalModel, t(testData))
  predClass <- rep("Intermediate", length(pred))
  predClass[pred < 3.5] <- "Low"
  predClass[pred >= 7.5] <- "High"
  
  obs <- Y_test$CAIDE
  obsClass <- rep("Intermediate", length(obs))
  obsClass[obs < 3.5] <- "Low"
  obsClass[obs >= 7.5] <- "High" 
  
  #accuracy[i] <- sum(predClass == obsClass)/length(predClass)
  accuracy[i] <- caret::confusionMatrix(factor(predClass, levels = c("Low", "Intermediate", "High")), 
                                        factor(obsClass, levels = c("Low", "Intermediate", "High")))$overall[2]
}


methods <- c("EN", "sPLSDA", "RF")
methodNames <- c("ElasticNet",  "sPLS-DA", "Random Forest")

for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_CAIDE1/Performance_CAIDE1_HighRisk_Cor_", methods[i],".RData"))
  load(paste0("CV_CAIDE1/Performance_CAIDE1_LowRisk_Cor_", methods[i],".RData"))
  
  # predicted in test set
  predClass <- rep("Intermediate",nrow(ObsPred_test_HighRisk))
  predClass[(ObsPred_test_HighRisk$pred == "High") &
                   (ObsPred_test_LowRisk$pred == "Intermediate_High")] <- "High"
  predClass[(ObsPred_test_LowRisk$pred == "Low") &
                   (ObsPred_test_HighRisk$pred == "Low_Intermediate")] <- "Low"
  
  # Observed in test set
  obs <- ObsPred_test_HighRisk$ObservedScore
  obsClass <- rep("Intermediate", length(obs))
  obsClass[obs < 3.5] <- "Low"
  obsClass[obs >= 7.5] <- "High"
  
  #accuracy[i+3] <- sum(predClass == obsClass)/length(predClass)
  accuracy[i+3] <- caret::confusionMatrix(factor(predClass, levels = c("Low", "Intermediate", "High")), 
                                          factor(obsClass, levels = c("Low", "Intermediate", "High")))$overall[2]
  
}


plotDF_CAIDE1 <- data.frame(Accuracy = accuracy,
                     Model = c(rep("Regression Model",3),
                               rep("Classification Model",3)),
                     Method = factor(rep(c("ElasticNet", "sPLS(-DA)", "Random Forest"),2),
                                     levels = c("ElasticNet", "sPLS(-DA)", "Random Forest")),
                     Score = rep("CAIDE1",6)
                     )

#==============================================================================#
# CAIDE2
#==============================================================================#
# Get test data
load("X/X_Cor_CAIDE2/X_test_Cor2.RData")
testData <- log2(X_test_Cor2/(1-X_test_Cor2))

accuracy <- rep(NA,6)
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

for (i in 1:length(methods)){
  # Load model
  load(paste0("CV_CAIDE2/CV_CAIDE2_Cor_", methods[i],".RData"))
  
  # make prediction
  pred <- predict(finalModel, t(testData))
  predClass <- rep("Intermediate", length(pred))
  predClass[pred < 4.5] <- "Low"
  predClass[pred >= 8.5] <- "High"
  
  obs <- Y_test$CAIDE
  obsClass <- rep("Intermediate", length(obs))
  obsClass[obs < 4.5] <- "Low"
  obsClass[obs >= 8.5] <- "High" 
  
  #accuracy[i] <- sum(predClass == obsClass)/length(predClass)
  accuracy[i] <- caret::confusionMatrix(factor(predClass, levels = c("Low", "Intermediate", "High")), 
                                        factor(obsClass, levels = c("Low", "Intermediate", "High")))$overall[2]
}


methods <- c("EN", "sPLSDA", "RF")
methodNames <- c("ElasticNet",  "sPLS-DA", "Random Forest")

for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_CAIDE2/Performance_CAIDE2_HighRisk_Cor_", methods[i],".RData"))
  load(paste0("CV_CAIDE2/Performance_CAIDE2_LowRisk_Cor_", methods[i],".RData"))
  
  # predicted in test set
  predClass <- rep("Intermediate",nrow(ObsPred_test_HighRisk))
  predClass[(ObsPred_test_HighRisk$pred == "High") &
              (ObsPred_test_LowRisk$pred == "Intermediate_High")] <- "High"
  predClass[(ObsPred_test_LowRisk$pred == "Low") &
              (ObsPred_test_HighRisk$pred == "Low_Intermediate")] <- "Low"
  
  # Observed in test set
  obs <- ObsPred_test_HighRisk$ObservedScore
  obsClass <- rep("Intermediate", length(obs))
  obsClass[obs < 4.5] <- "Low"
  obsClass[obs >= 8.5] <- "High"
  
  #accuracy[i+3] <- sum(predClass == obsClass)/length(predClass)
  accuracy[i+3] <- caret::confusionMatrix(factor(predClass, levels = c("Low", "Intermediate", "High")), 
                                          factor(obsClass, levels = c("Low", "Intermediate", "High")))$overall[2]
  
}


plotDF_CAIDE2 <- data.frame(Accuracy = accuracy,
                            Model = c(rep("Regression Model",3),
                                      rep("Classification Model",3)),
                            Method = factor(rep(c("ElasticNet", "sPLS(-DA)", "Random Forest"),2),
                                            levels = c("ElasticNet", "sPLS(-DA)", "Random Forest")),
                            Score = rep("CAIDE2",6)
)

#==============================================================================#
# LIBRA
#==============================================================================#
# Get test data
load("X/X_Cor_LIBRA/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))

accuracy <- rep(NA,6)
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

for (i in 1:length(methods)){
  # Load model
  load(paste0("CV_LIBRA/CV_LIBRA_Cor_", methods[i],".RData"))
  
  # make prediction
  pred <- predict(finalModel, t(testData))
  predClass <- rep("Intermediate", length(pred))
  predClass[pred < 0] <- "Low"
  predClass[pred > 2] <- "High"
  
  obs <- Y_test$LIBRA
  obsClass <- rep("Intermediate", length(obs))
  obsClass[obs < 0] <- "Low"
  obsClass[obs > 2] <- "High" 
  
  #accuracy[i] <- sum(predClass == obsClass)/length(predClass)
  accuracy[i] <- caret::confusionMatrix(factor(predClass, levels = c("Low", "Intermediate", "High")), 
                                        factor(obsClass, levels = c("Low", "Intermediate", "High")))$overall[2]
}


methods <- c("EN", "sPLSDA", "RF")
methodNames <- c("ElasticNet",  "sPLS-DA", "Random Forest")

for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_LIBRA/Performance_LIBRA_HighRisk_Cor_", methods[i],".RData"))
  load(paste0("CV_LIBRA/Performance_LIBRA_LowRisk_Cor_", methods[i],".RData"))
  
  # predicted in test set
  predClass <- rep("Intermediate",nrow(ObsPred_test_HighRisk))
  predClass[(ObsPred_test_HighRisk$pred == "High") &
              (ObsPred_test_LowRisk$pred == "Intermediate_High")] <- "High"
  predClass[(ObsPred_test_LowRisk$pred == "Low") &
              (ObsPred_test_HighRisk$pred == "Low_Intermediate")] <- "Low"
  
  # Observed in test set
  obs <- ObsPred_test_HighRisk$ObservedScore
  obsClass <- rep("Intermediate", length(obs))
  obsClass[obs < 0] <- "Low"
  obsClass[obs > 2] <- "High"
  
  #accuracy[i+3] <- sum(predClass == obsClass)/length(predClass)
  accuracy[i+3] <- caret::confusionMatrix(factor(predClass, levels = c("Low", "Intermediate", "High")), 
                                          factor(obsClass, levels = c("Low", "Intermediate", "High")))$overall[2]
  
}


plotDF_LIBRA <- data.frame(Accuracy = accuracy,
                            Model = c(rep("Regression Model",3),
                                      rep("Classification Model",3)),
                            Method = factor(rep(c("ElasticNet", "sPLS(-DA)", "Random Forest"),2),
                                            levels = c("ElasticNet", "sPLS(-DA)", "Random Forest")),
                            Score = rep("LIBRA",6)
)


plotDF <- rbind.data.frame(plotDF_CAIDE1, plotDF_CAIDE2, plotDF_LIBRA)
p <- ggplot(plotDF) +
  geom_bar(aes(x = Model, y = Accuracy, fill = Method), 
           stat = "identity", position = position_dodge(), color = "grey", linewidth = 1) +
  facet_grid(cols = vars(Score)) +
  scale_fill_manual(values = c("#EF3B2C","#FE9929", "#807DBA")) +
  xlab(NULL) +
  ylab("Cohen's kappa") +
  labs(fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(p, file = "Compare_CatCon.png", width = 8, height = 6)