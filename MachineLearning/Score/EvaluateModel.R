
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

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Meta data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

###############################################################################

# CAIDE1

###############################################################################

# Load data
load("CV_CAIDE1/CV_CAIDE1_Cor_EN.RData")
load("X/X_Cor_CAIDE1/X_test_Cor.RData")
testData <- log2(X_test_Cor/(1-X_test_Cor))

pred_test <- predict(finalModel, t(testData))
obs_test <- Y_test$CAIDE
rmse <- caret::RMSE(pred = pred_test,
            obs =obs_test)

plotDF <- data.frame(Predicted = pred_test,
                     Observed = obs_test)

ObsVsPred <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(x = Observed, y = Predicted, color = Observed-Predicted), size = 3, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("CAIDE1", subtitle = "ElasticNet-regularized linear regression model") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(ObsVsPred, file = "ObsVsPred_test_EN_CAIDE1.png", width = 8, height = 6)

#*****************************************************************************#
#* Which factors (R2)
#*****************************************************************************#

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/CAIDE.Rdata")
rownames(CAIDE) <- CAIDE$Basename
Y_test_CAIDE1 <- CAIDE[colnames(testData),]

# tested methods
methods <- c("EN")
methodNames <- c("ElasticNet")

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
  Y_factors <- Y_test_CAIDE1[,c(10:16)]
  for (factor in 1:ncol(Y_factors)){
    formula <- paste0("PredictedScore ~ ", colnames(Y_factors)[factor])
    
    # Fit model
    model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
    
    # Calculate R-squared
    Rsquared_factor[factor] = summary(model)$r.squared
  }
  
  factorNames <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                   "Total Cholesterol", "Physical Inactivity")
  plotDF <- data.frame(Name = factorNames, 
                      R2  = Rsquared_factor,
                       Method = rep(methodNames[i],length(factorNames)))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

# Make plot
p <- ggplot(plotDF_all)+
  geom_bar(aes(x = Name, y = R2, fill = Method),
           stat = "identity", position=position_dodge()) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression(R^2)) +
  ylim(c(0,1)) +
  ggtitle("CAIDE1") +
  labs(fill = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#FD841F", "#E14D2A", "#CD104D")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))


ggsave(p, file = "WhichFactors_CAIDE1.png", width = 8, height = 6)



###############################################################################

# CAIDE2

###############################################################################

# Load data
load("CV_CAIDE2/CV_CAIDE2_Cor_EN.RData")
load("X/X_Cor_CAIDE2/X_test_Cor2.RData")
testData <- log2(X_test_Cor2/(1-X_test_Cor2))

pred_test <- predict(finalModel, t(testData))
obs_test <- Y_test$CAIDE
rmse <- caret::RMSE(pred = pred_test,
                    obs =obs_test)

plotDF <- data.frame(Predicted = pred_test,
                     Observed = obs_test)

ObsVsPred <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(x = Observed, y = Predicted, color = Observed-Predicted), size = 3, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("CAIDE2", subtitle = "ElasticNet-regularized linear regression model") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(ObsVsPred, file = "ObsVsPred_test_EN_CAIDE2.png", width = 8, height = 6)

#*****************************************************************************#
#* Which factors (R2)
#*****************************************************************************#

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/CAIDE2.Rdata")
rownames(CAIDE2) <- CAIDE2$Basename
Y_test_CAIDE2 <- CAIDE2[colnames(testData),]

# tested methods
methods <- c("EN")
methodNames <- c("ElasticNet")

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
  Y_factors <- Y_test_CAIDE2[,c(10:17)]
  for (factor in 1:ncol(Y_factors)){
    formula <- paste0("PredictedScore ~ ", colnames(Y_factors)[factor])
    
    # Fit model
    model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
    
    # Calculate R-squared
    Rsquared_factor[factor] = summary(model)$r.squared
  }
  
  factorNames <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                   "Total Cholesterol", "Physical Inactivity", "APOE Status")
  plotDF <- data.frame(Name = factorNames, 
                       R2  = Rsquared_factor,
                       Method = rep(methodNames[i],length(factorNames)))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

# Make plot
p <- ggplot(plotDF_all)+
  geom_bar(aes(x = Name, y = R2, fill = Method),
           stat = "identity", position=position_dodge()) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression(R^2)) +
  ylim(c(0,1)) +
  ggtitle("CAIDE2") +
  labs(fill = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#FD841F", "#E14D2A", "#CD104D")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))


ggsave(p, file = "WhichFactors_CAIDE2.png", width = 8, height = 6)


#*****************************************************************************#
#* Which factors (cohenF2)
#*****************************************************************************#

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/CAIDE2.Rdata")
rownames(CAIDE2) <- CAIDE2$Basename
Y_test_CAIDE2 <- CAIDE2[colnames(testData),]

# tested methods
methods <- c("EN")
methodNames <- c("ElasticNet")

# For each method calculate cohen's f2 for each factor
plotDF_all <- NULL
for (i in 1:length(methods)){
  # Load data
  load(paste0("CV_CAIDE2/CV_CAIDE2_Cor_",methods[i],".RData"))
  
  # Make prediction
  pred_test <- predict(finalModel, t(testData))
  
  # Combine observed and predicted
  dataMatrix <- cbind.data.frame(Y_test_CAIDE2[,c(10:17)],pred_test)
  colnames(dataMatrix) <- c(colnames(Y_test_CAIDE2[,c(10:17)]), "PredictedScore")
  
  # Make formula
  formula <- paste0("PredictedScore ~ ", 
                    paste0("0 + ", paste(colnames(Y_test_CAIDE2[,c(10:17)]), collapse = " + ")))
  
  
  # Fit linear model
  model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
  
  # Get fitted and residual values
  fittedValues <- fitted(model)
  residualValues <- residuals(model)
  
  # Calculate R-squared
  sse <- sum(residualValues^2)
  ssr <- sum(fittedValues^2)
  Rsquared = ssr/(ssr + sse)
  
  # Get global effect size
  globalEffect <- (Rsquared)/(1-Rsquared)
  
  cohenF <- list()
  Y_factors <- Y_test_CAIDE2[,c(10:17)]
  for (factor in 1:ncol(Y_factors)){
    formula <- paste0("PredictedScore ~ ", 
                      paste0("0 + ", paste(colnames(Y_factors[,-factor]), collapse = " + ")))
    
    
    model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
    
    # Get fitted and residual values
    fittedValues <- fitted(model)
    residualValues <- residuals(model)
    
    # Calculate R-squared
    sse <- sum(residualValues^2)
    ssr <- sum(fittedValues^2)
    Rsquared_factor = ssr/(ssr + sse)
    
    # Calculate cohen's f2 statistic (local effect size)
    cohenF[[factor]] <- ((Rsquared - Rsquared_factor)/(1-Rsquared))
  }
  
  factorNames <- c("Age", "Sex", "Education", "Systolic Blood Pressure", "BMI",
                   "Total Cholesterol", "Physical Inactivity", "APOE")
  plotDF <- data.frame(Name = c(factorNames, "Global"), 
                       cohenF2  = c(unlist(cohenF), globalEffect),
                       Method = rep(methodNames[i],length(factorNames)+1))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

# Make plot
p <- ggplot(plotDF_all[plotDF$Name != "Global",])+
  geom_bar(aes(x = Name, y = cohenF2, fill = Method),
           stat = "identity", position=position_dodge()) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression("Cohen's " ~ f^2)) +
  ggtitle("CAIDE2") +
  labs(fill = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#FD841F", "#E14D2A", "#CD104D")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))


ggsave(p, file = "WhichFactors_CAIDE2.png", width = 8, height = 6)



###############################################################################

# LIBRA

###############################################################################

# load data
load("CV_LIBRA/CV_LIBRA_Cor_EN.RData")
load("X/X_Cor_LIBRA/X_test_CorL.RData")
testData <- log2(X_test_CorL/(1-X_test_CorL))


pred_test <- predict(finalModel, t(testData))
obs_test <- Y_test$LIBRA
rmse <- caret::RMSE(pred = pred_test,
                    obs =obs_test)

plotDF <- data.frame(Predicted = pred_test,
                     Observed = obs_test)

ObsVsPred <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(x = Observed, y = Predicted, color = Observed-Predicted), size = 3, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  #ggtitle("LIBRA", subtitle = "ElasticNet-regularized linear regression model") +
  ggtitle("LIBRA", subtitle = "sPLS regression model") +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(ObsVsPred, file = "ObsVsPred_test_sPLS_LIBRA.png", width = 8, height = 6)


#*****************************************************************************#
#* Which factors (R2)
#*****************************************************************************#

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/EPILIBRA.Rdata")
rownames(EPILIBRA) <- EPILIBRA$Basename
Y_test_LIBRA <- EPILIBRA[colnames(testData),]

# tested methods
methods <- c("EN", "sPLS")
methodNames <- c("ElasticNet", "sPLS")

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
  
  Rsquared_factor <- rep(NA, 10)
  Y_factors <- Y_test_LIBRA[,c(10:20)]
  for (factor in 1:ncol(Y_factors)){
    formula <- paste0("PredictedScore ~ ", colnames(Y_factors)[factor])
    
    # Fit model
    model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
    
    # Calculate R-squared
    Rsquared_factor[factor] = summary(model)$r.squared
  }
  
  factorNames <- c("Diet", "Physical Act.", "Smoking", "L-M Alcohol", "Obesity",
                   "Depression", "Diabetes", "Hypertension", "HDL", "Heart Dis.", "Kidney Dis.")
  plotDF <- data.frame(Name = factorNames, 
                       R2  = Rsquared_factor,
                       Method = rep(methodNames[i],length(factorNames)))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

# Make plot
p <- ggplot(plotDF_all)+
  geom_bar(aes(x = Name, y = R2, fill = Method),
           stat = "identity", position=position_dodge()) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression(R^2)) +
  ylim(c(0,1)) +
  ggtitle("LIBRA") +
  labs(fill = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#FD841F", "#E14D2A", "#CD104D")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))


ggsave(p, file = "WhichFactors_LIBRA.png", width = 8, height = 6)

#*****************************************************************************#
#* Which factors
#*****************************************************************************#

# Get factors of test set
load("E:/Thesis/EXTEND/Phenotypes/EPILIBRA.Rdata")
rownames(EPILIBRA) <- EPILIBRA$Basename
Y_test_LIBRA <- EPILIBRA[colnames(testData),]

# tested methods
methods <- c("EN", "sPLS")
methodNames <- c("ElasticNet", "sPLS")

# For each method calculate cohen's f2 for each factor
plotDF_all <- NULL
for (i in 1:length(methods)){
  # Load data
  load(paste0("CV_LIBRA/CV_LIBRA_Cor_",methods[i],".RData"))
  
  # Make prediction
  pred_test <- predict(finalModel, t(testData))
  
  # Combine observed and predicted
  dataMatrix <- cbind.data.frame(Y_test_LIBRA[,c(10:20)],pred_test)
  colnames(dataMatrix) <- c(colnames(Y_test_LIBRA[,c(10:20)]), "PredictedScore")
  
  # Make formula
  formula <- paste0("PredictedScore ~ ", 
                    paste0("0 + ", paste(colnames(Y_test_LIBRA[,c(10:20)]), collapse = " + ")))
  
  
  # Fit linear model
  model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
  
  # Get fitted and residual values
  fittedValues <- fitted(model)
  residualValues <- residuals(model)
  
  # Calculate R-squared
  sse <- sum(residualValues^2)
  ssr <- sum(fittedValues^2)
  Rsquared = ssr/(ssr + sse)
  
  # Get global effect size
  globalEffect <- (Rsquared)/(1-Rsquared)
  
  cohenF <- list()
  Y_factors <- Y_test_LIBRA[,c(10:20)]
  for (factor in 1:ncol(Y_factors)){
    formula <- paste0("PredictedScore ~ ", 
                      paste0("0 + ", paste(colnames(Y_factors[,-factor]), collapse = " + ")))
    
    
    model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
    
    # Get fitted and residual values
    fittedValues <- fitted(model)
    residualValues <- residuals(model)
    
    # Calculate R-squared
    sse <- sum(residualValues^2)
    ssr <- sum(fittedValues^2)
    Rsquared_factor = ssr/(ssr + sse)
    
    # Calculate cohen's f2 statistic (local effect size)
    cohenF[[factor]] <- ((Rsquared - Rsquared_factor)/(1-Rsquared))
  }
  
  factorNames <- c("Diet", "Physical Act.", "Smoking", "L-M Alcohol", "Obesity",
                   "Depression", "Diabetes", "Hypertension", "HDL", "Heart Dis.", "Kidney Dis.")
  plotDF <- data.frame(Name = c(factorNames, "Global"), 
                       cohenF2  = c(unlist(cohenF), globalEffect),
                       Method = rep(methodNames[i],length(factorNames)+1))
  
  plotDF_all <- rbind.data.frame(plotDF_all, plotDF)
}

# Make plot
p <- ggplot(plotDF_all[plotDF$Name != "Global",])+
  geom_bar(aes(x = Name, y = cohenF2, fill = Method),
           stat = "identity", position=position_dodge()) +
  coord_flip() +
  xlab(NULL) +
  ylab(expression("Cohen's " ~ f^2)) +
  ggtitle("CAIDE2") +
  labs(fill = NULL) +
  theme_minimal() +
  scale_fill_manual(values = c("#FD841F", "#E14D2A", "#CD104D")) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))


ggsave(p, file = "WhichFactors_LIBRA.png", width = 8, height = 6)

