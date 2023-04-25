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
files <- list.files('Y', ".RData")
for (f in files){
  load(paste0("Y/",f))
}

# Load data
load("E:/Thesis/EXTEND/Methylation/PerFactor/CombineFactors/predictedScore_factors_EXTEND.RData")

################################################################################

# CAIDE1

################################################################################

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/CAIDE.Rdata")
rownames(CAIDE) <- CAIDE$Basename
Y_test_CAIDE1 <- CAIDE[Y_test$Basename,]

# tested methods
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

# For each method calculate R2 for each factor
testData <- predictedScore_factors[Y_test$Basename,]
plotDF_all <- NULL
for (i in 1:length(methods)){
  # Load data
  load(paste0("E:/Thesis/EXTEND/Methylation/PerFactor/CombineFactors/Fit_CombineFactors_CAIDE1_",methods[i],".RData"))
  
  # Make prediction
  pred_test <- predict(fit, testData)
  
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


ggsave(p, file = "WhichFactors_CAIDE1_perfactor.png", width = 6, height = 4.5)


################################################################################

# CAIDE2

################################################################################

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/CAIDE2.Rdata")
rownames(CAIDE2) <- CAIDE2$Basename
Y_test_CAIDE2 <- CAIDE2[Y_test$Basename,]

# tested methods
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

# For each method calculate R2 for each factor
testData <- predictedScore_factors[Y_test$Basename,]
plotDF_all <- NULL
for (i in 1:length(methods)){
  # Load data
  load(paste0("E:/Thesis/EXTEND/Methylation/PerFactor/CombineFactors/Fit_CombineFactors_CAIDE2_",methods[i],".RData"))
  
  # Make prediction
  pred_test <- predict(fit, testData)
  
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


ggsave(p, file = "WhichFactors_CAIDE2_perfactor.png", width = 6, height = 4.5)

################################################################################

# LIBRA

################################################################################

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/EPILIBRA.Rdata")
rownames(EPILIBRA) <- EPILIBRA$Basename
Y_test_LIBRA <- EPILIBRA[Y_test$Basename,]

# tested methods
methods <- c("EN", "sPLS", "RF")
methodNames <- c("ElasticNet", "sPLS", "Random Forest")

# For each method calculate R2 for each factor
testData <- predictedScore_factors[Y_test$Basename,]
plotDF_all <- NULL
for (i in 1:length(methods)){
  # Load data
  load(paste0("E:/Thesis/EXTEND/Methylation/PerFactor/CombineFactors/Fit_CombineFactors_LIBRA_",methods[i],".RData"))
  
  # Make prediction
  pred_test <- predict(fit, testData)
  
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


ggsave(p, file = "WhichFactors_LIBRA_perfactor.png", width = 6, height = 4.5)