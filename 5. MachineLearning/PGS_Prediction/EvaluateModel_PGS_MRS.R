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
load("E:/Thesis/EXTEND/Methylation/PerFactor/CombineFactors/predictedScore_factors_EXTEND.RData")
load("E:/Thesis/EXTEND/df_list.RData")
PGS <- df_list$bayesr
colnames(PGS) <- paste0(colnames(PGS), "_PGS")
samples_test <- Y_test[Y_test$ID %in% rownames(PGS),]


###############################################################################

# All scores

###############################################################################

# Prepare test data
X_test <- cbind.data.frame(predictedScore_factors[samples_test$Basename,],
                           PGS[samples_test$ID,])
Y_test <- samples_test
all(rownames(X_test) == Y_test$Basename)

# Make predictions: MRSs + PGSs
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")
scores <- c("CAIDE1", "CAIDE2", "LIBRA")
plotDF_test <- NULL
for (s in 1:length(scores)){
  for (m in 1:length(methods)){
    # Load data
    if (file.exists(paste0("HybridPerFactor/Fit_CombineFactors_",scores[s],"_", methods[m],"_PGS.RData"))){
      load(paste0("HybridPerFactor/Fit_CombineFactors_",scores[s],"_", methods[m],"_PGS.RData"))
      
      if (scores[s] == "CAIDE1"){
        test <- predict(fit,X_test)
        perf_test <- R2(pred = test, obs = samples_test$CAIDE)
      }
      if (scores[s] == "CAIDE2"){
        test <- predict(fit,X_test)
        perf_test <- R2(pred = test, obs = samples_test$CAIDE2)
      }
      if (scores[s] == "LIBRA"){
        test <- predict(fit,X_test)
        perf_test <- R2(pred = test, obs = samples_test$LIBRA)
      }
      temp_test <- data.frame(R2 = perf_test,
                              Score = scores[s],
                              Method = methodNames[m])
      plotDF_test <- rbind.data.frame(plotDF_test, temp_test)
      
    }
    
  }
  
}

plotDF_test$Method <- factor(plotDF_test$Method,
                             levels = methodNames)

plotDF_test$Feature <- rep("MRS + PGS", nrow(plotDF_test))



# Prepare test data
X_test <- predictedScore_factors[samples_test$Basename,]
Y_test <- samples_test
all(rownames(X_test) == Y_test$Basename)

# Make predictions: MRSs only
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")
scores <- c("CAIDE1", "CAIDE2", "LIBRA")
plotDF1 <- NULL
plotDF_test1 <- NULL
for (s in 1:length(scores)){
  for (m in 1:length(methods)){
    # Load data
    if (file.exists(paste0("PerFactor/CombineFactors/Fit_CombineFactors_",scores[s],"_", methods[m],".RData"))){
      load(paste0("PerFactor/CombineFactors/Fit_CombineFactors_",scores[s],"_", methods[m],".RData"))
      
      if (scores[s] == "CAIDE1"){
        pred_test <- predict(fit, X_test)
        perf_test <- R2(pred = pred_test,
                        obs = Y_test$CAIDE)
      }
      if (scores[s] == "CAIDE2"){
        pred_test <- predict(fit, X_test)
        perf_test <- R2(pred = pred_test,
                        obs = Y_test$CAIDE2)
      }
      if (scores[s] == "LIBRA"){
        pred_test <- predict(fit, X_test)
        perf_test <- R2(pred = pred_test,
                        obs = Y_test$LIBRA)
      }
      temp_test <- data.frame(R2 = perf_test,
                              Score = scores[s],
                              Method = methodNames[m])
      plotDF_test1 <- rbind.data.frame(plotDF_test1, temp_test)
      
    }
    
  }
  
}

plotDF_test1$Method <- factor(plotDF_test1$Method,
                              levels = methodNames)

plotDF_test1$Feature <- rep("MRS", nrow(plotDF_test1))


# Prepare data for plotting
plotDF_test_all <- rbind.data.frame(plotDF_test, plotDF_test1)
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

# Set colors
colors <- c(RColorBrewer::brewer.pal(n = 8, name = "Reds")[6],
            RColorBrewer::brewer.pal(n = 8, name = "Oranges")[6],
            RColorBrewer::brewer.pal(n = 8, name = "PuRd")[6])

# Make plot
p <- ggplot(testDF) +
  geom_bar(aes(y = R2, x = Method, fill = Score, alpha = Feature), linewidth = 0.7, 
           position = position_dodge(), stat = "identity", color = "black") +
  facet_grid(cols = vars(Score), scale = "free", space = "free") +
  guides(fill = "none") +
  ylab(expression(R^2)) +
  xlab(NULL)+
  scale_alpha_manual(values = c(0.5,1)) +
  scale_fill_manual(values = c("#FB6A4A","#FEC44F", "#BCBDDC")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# Save plot
ggsave(p,file = "Evaluate_MRSPGS_test1.png", width = 8, height = 5)


# Make alternative plot
p <- ggplot(testDF) +
  #geom_segment(aes(y = 0, yend = R2, x = Method, xend = Method), linewidth = 1, color = "black") +
  geom_line(aes(y = R2, x = Method, group = Group, color = Color), linewidth = 1.5) +
  geom_point(aes(y = R2,x = Method, shape = Feature, size = Feature, color = Color1)) +
  facet_grid(cols = vars(Score), scale = "free", space = "free") +
  guides(color = "none") +
  xlab(expression(R^2)) +
  ylab(NULL)+
  scale_size_manual(values = c(3,3)) +
  scale_shape_manual(values = c(15,16)) +
  scale_color_manual(values = c("red","black", "red", "blue","black", "blue")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# Save alternative plot
ggsave(p,file = "Evaluate_MRSPGS_test.png", width = 8, height = 5)


