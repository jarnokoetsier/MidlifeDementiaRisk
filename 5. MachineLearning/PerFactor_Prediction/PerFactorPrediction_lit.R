
# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

################################################################################

# Model Training

################################################################################

# Load data
load("Data/X_test.RData")
load("Data/X_nonTest.RData")
load("PerFactor/Y_test_factors.RData")
load("PerFactor/Y_nonTest_factors.RData")
load("PerFactor/selectedProbes_lit.RData")

# Test if samples are in correct order
Y_train <- Y_nonTest[,c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", 
                        "Physical", "HeartDisease", "Education", "Depression",
                        "SysBP", "Diet")]
X_train <- log2(X_nonTest/(1-X_nonTest))

Y_test <- Y_test[,c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", 
                    "Physical", "HeartDisease", "Education", "Depression",
                    "SysBP", "Diet")]
X_test <- log2(X_test/(1-X_test))

all(colnames(X_train) == rownames(Y_train))
all(colnames(X_test) == rownames(Y_test))

# Set number of folds and repeats
nfold = 5
nrep = 5

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Machine learning method
MLmethod = "glmnet"

outputList <- list()
for (f in 1:ncol(Y_train)){
  
  
  # Performance metric
  performance_metric = "ROC"
  
  # factor
  factor1 <- Y_train[!is.na(Y_train[,f]),f]
  
  # dataMatrix
  dataMatrix <- X_train[rownames(X_train) %in% selectedProbes_lit[[f]],!is.na(Y_train[,f])]
  
  # Create index
  set.seed(123)
  CVindex <- NULL
  for (r in 1:5){
    temp <- createFolds(1:length(factor1),5, returnTrain = TRUE)
    names(temp) <- paste0(names(temp),".Rep",r)
    CVindex <- c(CVindex, temp)
  }
  rm(temp)
  
  # Model training
  trainResults <- list()
  for(i in 1:length(CVindex)) {
    
    # Sample indices for repeated CV
    index <- list(CVindex[[i]])
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = FALSE,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               index = index)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix),
                 y = factor(factor1, levels = c("No", "Yes")),
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)
    
    
    trainResults[[i]] <- fit$results
    
  }
  
  # Save output
  outputList[[f]] <- trainResults
  
  
}

names(outputList) <- colnames(Y_train)
save(outputList, file = "PerFactor/OutputList_factors_lit.RData")


################################################################################

# Model Performance

################################################################################
load("PerFactor/OutputList_factors_lit.RData")

PerformanceList <- list()
for (i in 1:length(outputList)){
  trainResults <- outputList[[i]]
  perf <- matrix(NA, nrow = 1000, ncol = 25)
  for (j in 1:length(trainResults)){
    perf[,j] <- trainResults[[j]]$ROC
  }
  optPar <- which.max(rowMeans(perf))
  optPerf <- perf[optPar,]
  optAlpha <- trainResults[[1]]$alpha[optPar]
  optLambda <- trainResults[[1]]$lambda[optPar]
  
  PerformanceList[[i]] <- list(optPerf, optAlpha, optLambda)
}
names(PerformanceList) <- names(outputList)
save(PerformanceList, file = "PerFactor/PerformanceList_factors_Lit.RData")

################################################################################

# Obs vs pred in CV

################################################################################

# Load data
load("Data/X_test.RData")
load("Data/X_nonTest.RData")
load("PerFactor/Y_test_factors.RData")
load("PerFactor/Y_nonTest_factors.RData")
load("PerFactor/selectedProbes_lit.RData")

# Test if samples are in correct order
Y_train <- Y_nonTest[,c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", 
                        "Physical", "HeartDisease", "Education", "Depression")]
X_train <- log2(X_nonTest/(1-X_nonTest))

Y_test <- Y_test[,c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", 
                    "Physical", "HeartDisease", "Education", "Depression")]
X_test <- log2(X_test/(1-X_test))

all(colnames(X_train) == rownames(Y_train))
all(colnames(X_test) == rownames(Y_test))

load("PerFactor/PerformanceList_factors_Cor.RData")

ObsPred_all <- list()
for (f in 1:length(PerformanceList)){
  
  # Set grid for lambda
  lambdaCV <- PerformanceList[[f]][[3]]
  
  # Set grid for alpha
  alphaCV <- PerformanceList[[f]][[2]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(alphaCV, lambdaCV)
  colnames(parameterGrid) <- c(".alpha", ".lambda")
  
  # ML method
  MLmethod = "glmnet"
  
  # Performance metric
  performance_metric = "ROC"
  
  # factor
  factor1 <- Y_train[!is.na(Y_train[,f]),f]
  
  # dataMatrix
  dataMatrix <- X_train[rownames(X_train) %in% selectedProbes_lit[[f]],!is.na(Y_train[,f])]
  
  
  # Create index
  set.seed(123)
  CVindex <- NULL
  for (r in 1:5){
    temp <- createFolds(1:length(factor1),5, returnTrain = TRUE)
    names(temp) <- paste0(names(temp),".Rep",r)
    CVindex <- c(CVindex, temp)
  }
  rm(temp)
  
  # Model training
  ObsPred_CV <- NULL
  for(i in 1:length(CVindex)) {
    
    # Sample indices for repeated CV
    index <- list(CVindex[[i]])
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = TRUE,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               index = index)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix),
                 y = factor(factor1, levels = c("No", "Yes")),
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)
    
    ObsPred_CV <- rbind.data.frame(ObsPred_CV, fit$pred[,1:3])
  }
  ObsPred_all[[f]] <- ObsPred_CV
}
names(ObsPred_all) <- colnames(Y_train)
save(ObsPred_all, file = "PerFactor/ObsPred_lit.RData")


################################################################################

# finalOutput CV

################################################################################

# construct ROC
finalOutput <- list()
for (f in 1:length(ObsPred_all)){
  ObsPred_CV <- ObsPred_all[[f]]
  roc_list_CV <- roc(response = factor(ObsPred_CV$obs, levels = c("No", "Yes")), 
                     predictor = ObsPred_CV$No)
  
  # Get sensitivies and specificities
  plotDF <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                       Specificity = roc_list_CV$specificities)
  
  # Get AUC
  auc_CV <- auc(roc_list_CV)
  
  # Get predictions
  Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
  threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
  ObsPred_CV$PredictedClass <- ifelse(ObsPred_CV$No > threshold, "No", "Yes")
  
  # Save
  finalOutput[[f]] <- list(AUCplot = plotDF,
                           AUC = as.numeric(auc_CV),
                           threshold = threshold,
                           ObsPred_CV = ObsPred_CV)
}
names(finalOutput) <- names(ObsPred_all)
save(finalOutput, file = "PerFactor/finalOutput_lit_CV.RData")

factorNames <- c("BMI", "Diabetes", "L-M Alcohol", "HDL Chol", "Total Chol.", 
                 "Physical Act.", "Heart Disease", "Education", "Depression",
                 "Systolic Blood Pressure", "Dietary Intake")
p <- ggplot()
for (f in 1:length(finalOutput)){
  plotData <- finalOutput[[f]]$AUCplot
  plotData$Factor <- rep(factorNames[f],nrow(plotData))
  
  p <- p + 
    geom_path(data = plotData,  aes(y = Sensitivity, x = 1- Specificity, color = Factor),
              linewidth = 1)
}

plotCV <- p +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  labs(color = "Risk Factor") +
  theme_classic() +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 8, name = "Set2"),
                                RColorBrewer::brewer.pal(n = 8, name = "Dark2")))

ggsave(plotCV, file = "PerFactor/CVperformance_EN_perFactor.png", width = 8, height = 6)




################################################################################

# Obs vs pred test set

################################################################################


# Load data
load("Data/X_test.RData")
load("Data/X_nonTest.RData")
load("PerFactor/Y_test_factors.RData")
load("PerFactor/Y_nonTest_factors.RData")
load("PerFactor/selectedProbes_lit.RData")

# Test if samples are in correct order
Y_train <- Y_nonTest[,c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", 
                        "Physical", "HeartDisease", "Education", "Depression")]
X_train <- log2(X_nonTest/(1-X_nonTest))

Y_test <- Y_test[,c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", 
                    "Physical", "HeartDisease", "Education", "Depression")]
X_test <- log2(X_test/(1-X_test))

all(colnames(X_train) == rownames(Y_train))
all(colnames(X_test) == rownames(Y_test))

load("PerFactor/PerformanceList_factors_Cor.RData")


finalOutput_test <- list()
for (f in 1:length(PerformanceList)){
  
  # Set grid for lambda
  lambdaCV <- PerformanceList[[f]][[3]]
  
  # Set grid for alpha
  alphaCV <- PerformanceList[[f]][[2]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(alphaCV, lambdaCV)
  colnames(parameterGrid) <- c(".alpha", ".lambda")
  
  # ML method
  MLmethod = "glmnet"
  
  # Performance metric
  performance_metric = "ROC"
  
  # factor
  factor1 <- Y_train[!is.na(Y_train[,f]),f]
  factor_cor <- ifelse(factor1 == "Yes",1,0)
  
  # dataMatrix
  dataMatrix <- X_train[rownames(X_train) %in% selectedProbes_lit[[f]],!is.na(Y_train[,f])]
  
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(method = "none", 
                             classProbs = TRUE)
  
  # Actual training
  set.seed(123)
  finalModel <- train(x = t(dataMatrix),
                      y = factor(factor1, levels = c("No", "Yes")),
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = TRUE)
  
  # Get test predictions
  pred <- predict(finalModel, t(X_test[rownames(X_test) %in% selectedProbes_lit[[f]],]), type = "prob")
  obs <- Y_test[,f]
  threshold <- finalOutput[[f]]$threshold
  predClass <- ifelse(pred$No > threshold, "No", "Yes")
  
  ObsPred_test <- data.frame(pred = pred$No,
                             predClass = predClass,
                             obs = obs)
  
  # Get AUC
  roc_list_test <- roc(response = factor(ObsPred_test$obs, levels = c("No", "Yes")), 
                       predictor = ObsPred_test$pred)
  
  # Get sensitivies and specificities
  plotDF <- data.frame(Sensitivity = roc_list_test$sensitivities,
                       Specificity = roc_list_test$specificities)
  
  # Get AUC
  auc_test <- auc(roc_list_test)
  
  # Save
  finalOutput_test[[f]] <- list(AUCplot = plotDF,
                                AUC = as.numeric(auc_test),
                                threshold = threshold,
                                ObsPred_test = ObsPred_test)
  
}
names(finalOutput_test) <- colnames(Y_train)
save(finalOutput_test, file = "PerFactor/finalOutput_lit_test.RData")

# Performance on test set


factorNames <- c("BMI", "Diabetes", "L-M Alcohol", "HDL Chol", "Total Chol.", 
                 "Physical Act.", "Heart Disease", "Education", "Depression",
                 "Systolic Blood Pressure", "Dietary Intake")
p <- ggplot()
for (f in 1:length(finalOutput_test)){
  plotData <- finalOutput_test[[f]]$AUCplot
  plotData$Factor <- rep(factorNames[f],nrow(plotData))
  
  p <- p + 
    geom_path(data = plotData,  aes(y = Sensitivity, x = 1- Specificity, color = Factor),
              linewidth = 1)
}

plotTest <- p +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  labs(color = "Risk Factor") +
  theme_classic() +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 8, name = "Set2"),
                                RColorBrewer::brewer.pal(n = 8, name = "Dark2")))

ggsave(plotTest, file = "PerFactor/Testperformance_EN_perFactor.png", width = 8, height = 6)

