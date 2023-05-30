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

# Get test data
load("X/X_Cor_CAIDE1/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))

predNames <- c("CAIDE1 Score: ElasticNet",
               "CAIDE1 Score: sPLS",
               "CAIDE1 Score: Random Forest",
               "CAIDE1 Low Risk: ElasticNet",
               "CAIDE1 Low Risk: sPLS-DA",
               "CAIDE1 Low Risk: Random Forest",
               "CAIDE1 High Risk: ElasticNet",
               "CAIDE1 High Risk: sPLS-DA",
               "CAIDE1 High Risk: Random Forest")

predMatrix <- matrix(NA,nrow(Y_test), length(predNames))
colnames(predMatrix) <- predNames
rownames(predMatrix) <- Y_test$Basename

methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

for (i in 1:length(methods)){
  # Load model
  load(paste0("CV_CAIDE1/CV_CAIDE1_Cor_", methods[i],".RData"))
  
  # make prediction
  predMatrix[,i] <- predict(finalModel, t(testData))

}

methods <- c("EN", "sPLSDA", "RF")
methodNames <- c("ElasticNet",  "sPLS-DA", "Random Forest")

for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_CAIDE1/Performance_CAIDE1_LowRisk_Cor_", methods[i],".RData"))
  
  # get prediction
  predMatrix[,i+3] <- 1-ObsPred_test_LowRisk$p

  
}

methods <- c("EN", "sPLSDA", "RF")
methodNames <- c("ElasticNet",  "sPLS-DA", "Random Forest")

for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_CAIDE1/Performance_CAIDE1_HighRisk_Cor_", methods[i],".RData"))
  
  # get prediction
  predMatrix[,i+6] <- 1-ObsPred_test_HighRisk$p
  
}

predMatrix_scaled <- t((t(predMatrix) - rowMeans(t(predMatrix)))/apply(t(predMatrix),1,sd))

colMeans(predMatrix_scaled)
apply(predMatrix_scaled,2,sd)


plotPred <- gather(as.data.frame(predMatrix_scaled))
plotPred$Basename <- rep(rownames(predMatrix_scaled), ncol(predMatrix_scaled))



plotPred$Score <- rep(NA, nrow(plotPred))
plotPred$Score[str_detect(plotPred$key, "Score")] <- "CAIDE1 Score"
plotPred$Score[str_detect(plotPred$key, "Low Risk")] <- "Low CAIDE1"
plotPred$Score[str_detect(plotPred$key, "High Risk")] <- "High CAIDE1"

plotPred$Model <- rep("sPLS", nrow(plotPred))
plotPred$Model[str_detect(plotPred$key, "sPLS-DA")] <- "sPLS-DA"
plotPred$Model[str_detect(plotPred$key, "Random Forest")] <- "Random Forest"
plotPred$Model[str_detect(plotPred$key, "ElasticNet")] <- "ElasticNet"


plotObs <- Y_test[,c("Basename", "CAIDE")]
plotPred <- inner_join(plotPred, plotObs, by = c("Basename" = "Basename"))
plotPred$Cat <- rep("Intermediate", nrow(plotPred)) 
plotPred$Cat[plotPred$CAIDE < 4] <- "Low"
plotPred$Cat[plotPred$CAIDE > 7] <- "High"


sampleOrder <- Y_test$CAIDE
names(sampleOrder) <- Y_test$Basename
sampleOrder <- sort(sampleOrder)

plotPred$Basename <- factor(plotPred$Basename,
                            levels = names(sampleOrder))

main <- ggplot(plotPred) +
  geom_tile(aes(x = Basename, y = Model, fill = value),
            width = 1, height = 1.2) +
  facet_grid(rows = vars(Score), 
             scales = "free", space = "free") +
  scale_fill_viridis_c() +
  xlab("Samples") +
  ylab(NULL) +
  labs(fill = "Unit Scaled\nPredicted Value") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "right",
        axis.text.y = element_text(),
        strip.background = element_rect(fill = "grey", linewidth = 0),
        strip.text.y = element_text(angle = 270, color = "black"),
        panel.spacing.x = unit(0.2, "lines"))


plotObs$Basename <- factor(plotObs$Basename,
                           levels = names(sampleOrder))

top <- ggplot(plotObs) +
  geom_bar(aes(x = Basename, y = CAIDE, fill = CAIDE), width = 1, 
           stat = "identity") +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_gradient(low = "#FEE0D2", high = "#CB181D")



p <- top + main +
  plot_layout(nrow = 2, ncol = 1, height = c(1,5))

ggsave(p, file = "CompareModels_CAIDE1.png", width = 8, height = 6)






corMatrix <- correlate(predMatrix, method = "spearman")







methods <- c("EN", "sPLSDA", "RF")
methodNames <- c("ElasticNet",  "sPLS-DA", "Random Forest")

for (i in 1:length(methods)){
  
  # load data
  load(paste0("CV_CAIDE1/Performance_CAIDE1_HighRisk_Cor_", methods[i],".RData"))
  load(paste0("CV_CAIDE1/Performance_CAIDE1_LowRisk_Cor_", methods[i],".RData"))
  
  # Obs vs predicted in test set
  temp <- data.frame(Predicted = rep(2,nrow(ObsPred_test_HighRisk)),
                     Observed = ObsPred_test_HighRisk$ObservedScore)
  temp$Predicted[(ObsPred_test_HighRisk$pred == "High") &
               (ObsPred_test_LowRisk$pred == "Intermediate_High")] <- 3
  temp$Predicted[(ObsPred_test_LowRisk$pred == "Low") &
               (ObsPred_test_HighRisk$pred == "Low_Intermediate")] <- 1
  temp$Method <- rep(methodNames[i], nrow(temp))

  predMatrix[,i+3] <- temp$Predicted
  
}

corMatrix <- correlate(predMatrix, method = "spearman")


