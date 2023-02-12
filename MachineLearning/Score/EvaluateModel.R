
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

methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")
scores <- c("CAIDE1", "CAIDE2", "LIBRA")
plotDF <- NULL
plotDF_test <- NULL
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

p <- ggplot(plotDF) +
  geom_boxplot(aes(x = Score, y = R2, fill = Method), 
               alpha = 1, outlier.shape = NA) +
  geom_point(aes(x = Score, y = R2, color = Method), 
             position=position_jitterdodge(jitter.width = 0.1), 
             size = 2, alpha = 0.8) +
  geom_point(data = plotDF_test, aes(x = Score, y = R2, color = Method), 
            fill = "black", position=position_jitterdodge(jitter.width = 0), 
            size = 4, shape = 23, alpha = 0.7) +
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

ggsave(p,file = "RiskScores_test.png", width = 8, height = 6)

###############################################################################

# CAIDE1

###############################################################################

# Get test data
load("X/X_Cor_CAIDE1/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))

methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")
plotDF <-  NULL
for (i in 1:length(methods)){
  # Load model
  load(paste0("CV_CAIDE1/CV_CAIDE1_Cor_", methods[i],".RData"))

  # make prediction
  pred_test <- predict(finalModel, t(testData))
  
  # Get observed value
  obs_test <- Y_test$CAIDE

  # Combine into data frame
  temp <- data.frame(Predicted = pred_test,
                     Observed = obs_test,
                     Method = rep(methodNames[i], length(obs_test)))
  
  plotDF <- rbind.data.frame(plotDF, temp)
}
plotDF$Method <- factor(plotDF$Method,
                        levels = methodNames)

p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(x = Observed, y = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  facet_grid(rows = vars(Method)) +
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

ggsave(g, file = "ObsVsPred_test_EN_CAIDE1.png", width = 8, height = 8)

#*****************************************************************************#
#* Which factors (R2)
#*****************************************************************************#

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/CAIDE.Rdata")
rownames(CAIDE) <- CAIDE$Basename
Y_test_CAIDE1 <- CAIDE[colnames(testData),]

# tested methods
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

# For each method calculate R2 for each factor
plotDF_all <- NULL
for (i in 1:length(methods)){
  # Load data
  load(paste0("CV_CAIDE1/CV_CAIDE1_Cor_",methods[i],".RData"))
  
  # Make prediction
  pred_test <- predict(finalModel, t(testData))
  
  # Combine observed and predicted
  dataMatrix <- cbind.data.frame(Y_test_CAIDE1[,c(10:16)],pred_test)
  colnames(dataMatrix) <- c(colnames(Y_test_CAIDE1[,c(10:16)]), "PredictedScore")

  Rsquared_factor <- rep(NA, 7)
  Sig_factor <- rep(NA, 7)
  Y_factors <- Y_test_CAIDE1[,c(10:16)]
  for (factor in 1:ncol(Y_factors)){
    formula <- paste0("PredictedScore ~ ", colnames(Y_factors)[factor])
    
    # Fit model
    model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
    
    # Calculate R-squared
    Rsquared_factor[factor] <-  summary(model)$r.squared
    Sig_factor[factor] <- summary(model)$coefficients[2,4]
  }
  
  factorNames <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                   "Total Cholesterol", "Physical Inactivity")
  plotDF <- data.frame(Name = factorNames, 
                       R2  = Rsquared_factor,
                       Sig = Sig_factor,
                       Method = rep(methodNames[i],length(factorNames)))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

plotDF_all$Method <- factor(plotDF_all$Method,
                        levels = methodNames)


plotDF_all$X <- plotDF_all$R2 + ifelse(plotDF_all$Sig < 0.05, 0.01,NA)
# Make plot
p <- ggplot(plotDF_all)+
  geom_bar(aes(x = Name, y = R2, fill = Method),
           stat = "identity", position=position_dodge(), color = "black") +
  geom_point(aes(x = Name, y = X, color = Method),
             position=position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.9),
             shape = 18) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression(R^2)) +
  ylim(c(0,1)) +
  ggtitle("CAIDE1") +
  labs(fill = NULL) +
  guides(color = "none") +
  theme_minimal() +
  scale_fill_manual(values = c("#EF3B2C","#FE9929", "#807DBA")) +
  scale_color_manual(values = c("black","black", "black")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))


ggsave(p, file = "WhichFactors_CAIDE1.png", width = 8, height = 6)



###############################################################################

# CAIDE2

###############################################################################

# Get test data
load("X/X_Cor_CAIDE2/X_test_Cor2.RData")
testData <- log2(X_test_Cor2/(1-X_test_Cor2))

methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")
plotDF <-  NULL
for (i in 1:length(methods)){
  # Load model
  load(paste0("CV_CAIDE2/CV_CAIDE2_Cor_", methods[i],".RData"))
  
  # make prediction
  pred_test <- predict(finalModel, t(testData))
  
  # Get observed value
  obs_test <- Y_test$CAIDE
  
  # Combine into data frame
  temp <- data.frame(Predicted = pred_test,
                     Observed = obs_test,
                     Method = rep(methodNames[i], length(obs_test)))
  
  plotDF <- rbind.data.frame(plotDF, temp)
}
plotDF$Method <- factor(plotDF$Method,
                        levels = methodNames)

p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(x = Observed, y = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  facet_grid(rows = vars(Method)) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("CAIDE2") +
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

ggsave(g, file = "ObsVsPred_test_EN_CAIDE2.png", width = 8, height = 8)

#*****************************************************************************#
#* Which factors (R2)
#*****************************************************************************#

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/CAIDE2.Rdata")
rownames(CAIDE2) <- CAIDE2$Basename
Y_test_CAIDE2 <- CAIDE2[colnames(testData),]

# tested methods
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

# For each method calculate R2 for each factor
plotDF_all <- NULL
for (i in 1:length(methods)){
  # Load data
  load(paste0("CV_CAIDE2/CV_CAIDE2_Cor_",methods[i],".RData"))
  
  # Make prediction
  pred_test <- predict(finalModel, t(testData))
  
  # Combine observed and predicted
  dataMatrix <- cbind.data.frame(Y_test_CAIDE2[,c(10:17)],pred_test)
  colnames(dataMatrix) <- c(colnames(Y_test_CAIDE2[,c(10:17)]), "PredictedScore")
  
  Rsquared_factor <- rep(NA, 8)
  Sig_factor <- rep(NA, 8)
  Y_factors <- Y_test_CAIDE2[,c(10:17)]
  for (factor in 1:ncol(Y_factors)){
    formula <- paste0("PredictedScore ~ ", colnames(Y_factors)[factor])
    
    # Fit model
    model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
    
    # Calculate R-squared
    Rsquared_factor[factor] <-  summary(model)$r.squared
    Sig_factor[factor] <- summary(model)$coefficients[2,4]
  }
  
  factorNames <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                   "Total Cholesterol", "Physical Inactivity", "APOE")
  plotDF <- data.frame(Name = factorNames, 
                       R2  = Rsquared_factor,
                       Sig = Sig_factor,
                       Method = rep(methodNames[i],length(factorNames)))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

plotDF_all$Method <- factor(plotDF_all$Method,
                            levels = methodNames)


plotDF_all$X <- plotDF_all$R2 + ifelse(plotDF_all$Sig < 0.05, 0.01,NA)
# Make plot
p <- ggplot(plotDF_all)+
  geom_bar(aes(x = Name, y = R2, fill = Method),
           stat = "identity", position=position_dodge(), color = "black") +
  geom_point(aes(x = Name, y = X, color = Method),
             position=position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.9),
             shape = 18) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression(R^2)) +
  ylim(c(0,1)) +
  ggtitle("CAIDE2") +
  labs(fill = NULL) +
  guides(color = "none") +
  theme_minimal() +
  scale_fill_manual(values = c("#EF3B2C","#FE9929", "#807DBA")) +
  scale_color_manual(values = c("black","black", "black")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))


ggsave(p, file = "WhichFactors_CAIDE2.png", width = 8, height = 6)





###############################################################################

# LIBRA

###############################################################################

# Get test data
load("X/X_Cor_LIBRA/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))

methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")
plotDF <-  NULL
for (i in 1:length(methods)){
  # Load model
  load(paste0("CV_LIBRA/CV_LIBRA_Cor_", methods[i],".RData"))
  
  # make prediction
  pred_test <- predict(finalModel, t(testData))
  
  # Get observed value
  obs_test <- Y_test$LIBRA
  
  # Combine into data frame
  temp <- data.frame(Predicted = pred_test,
                     Observed = obs_test,
                     Method = rep(methodNames[i], length(obs_test)))
  
  plotDF <- rbind.data.frame(plotDF, temp)
}
plotDF$Method <- factor(plotDF$Method,
                        levels = methodNames)

p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(x = Observed, y = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  facet_grid(rows = vars(Method)) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("LIBRA") +
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

ggsave(g, file = "ObsVsPred_test_LIBRA.png", width = 8, height = 8)

#*****************************************************************************#
#* Which factors (R2)
#*****************************************************************************#

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/EPILIBRA.Rdata")
rownames(EPILIBRA) <- EPILIBRA$Basename
Y_test_LIBRA <- EPILIBRA[colnames(testData),]

# tested methods
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

# For each method calculate R2 for each factor
plotDF_all <- NULL
for (i in 1:length(methods)){
  # Load data
  load(paste0("CV_LIBRA/CV_LIBRA_Cor_",methods[i],".RData"))
  
  # Make prediction
  pred_test <- predict(finalModel, t(testData))
  
  # Combine observed and predicted
  dataMatrix <- cbind.data.frame(Y_test_LIBRA[,c(10:20)],pred_test)
  colnames(dataMatrix) <- c(colnames(Y_test_LIBRA[,c(10:20)]), "PredictedScore")
  
  Rsquared_factor <- rep(NA, 11)
  Sig_factor <- rep(NA, 11)
  Y_factors <- Y_test_LIBRA[,c(10:20)]
  for (factor in 1:ncol(Y_factors)){
    formula <- paste0("PredictedScore ~ ", colnames(Y_factors)[factor])
    
    # Fit model
    model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
    
    # Calculate R-squared
    Rsquared_factor[factor] <-  summary(model)$r.squared
    Sig_factor[factor] <- summary(model)$coefficients[2,4]
  }
  
  factorNames <- c("Healthy Diet", "Physical Inactivity", "Smoking", "L-M Alcohol Intake",
                   "BMI", "Depression", "Type 2 Diabetes","Systolic Blood Pressure", 
                   "HDL Cholesterol", "Heart Disease", "Kidney Disease")
  plotDF <- data.frame(Name = factorNames, 
                       R2  = Rsquared_factor,
                       Sig = Sig_factor,
                       Method = rep(methodNames[i],length(factorNames)))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

plotDF_all$Method <- factor(plotDF_all$Method,
                            levels = methodNames)


plotDF_all$X <- plotDF_all$R2 + ifelse(plotDF_all$Sig < 0.05, 0.01,NA)
# Make plot
p <- ggplot(plotDF_all)+
  geom_bar(aes(x = Name, y = R2, fill = Method),
           stat = "identity", position=position_dodge(), color = "black") +
  geom_point(aes(x = Name, y = X, color = Method),
             position=position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.9),
             shape = 18) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression(R^2)) +
  ylim(c(0,1)) +
  ggtitle("LIBRA") +
  labs(fill = NULL) +
  guides(color = "none") +
  theme_minimal() +
  scale_fill_manual(values = c("#EF3B2C","#FE9929", "#807DBA")) +
  scale_color_manual(values = c("black","black", "black")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))


ggsave(p, file = "WhichFactors_LIBRA.png", width = 8, height = 6)

