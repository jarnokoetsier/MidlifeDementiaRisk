# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(ranger)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)

# Load data
load("Data/metaData_ageFil.RData")
load("Data/X_nonTest.RData")
load("Data/X_test.RData")


################################################################################

# Model Training

################################################################################

# Load data
load("Data/X_test.RData")
load("Data/X_nonTest.RData")
load("PerFactor/Y_test_factors.RData")
load("PerFactor/Y_nonTest_factors.RData")
load("PerFactor/selectedProbes.RData")
selectedProbes <- selectedProbes[-13]

# Test if samples are in correct order
Y_train <- Y_nonTest[,-13] # ignore kidney disease
Y_test <- Y_test[,-13] # ignore kidney disease
X_train <- log2(X_nonTest/(1-X_nonTest))
X_test <- log2(X_test/(1-X_test))
all(colnames(X_train) == rownames(Y_train))
all(colnames(X_test) == rownames(Y_test))


# Set number of folds and repeats
nfold = 5
nrep = 5

# Number of randomly selected predictors
mtry_CV <- c(100, 500,1000,1500,2000,3000,4000)

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = c(3,5,10,15,20)

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
MLmethod = "ranger"


outputList <- list()
for (f in 1:ncol(Y_train)){
  
  
  # Performance metric
  performance_metric = "ROC"
  
  # factor
  factor <- Y_train[!is.na(Y_train[,f]),f]
  
  # dataMatrix
  dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
  
  # Create index
  set.seed(123)
  CVindex <- NULL
  for (r in 1:5){
    temp <- createFolds(1:length(factor),5, returnTrain = TRUE)
    names(temp) <- paste0(names(temp),".Rep",r)
    CVindex <- c(CVindex, temp)
  }
  rm(temp)
  
  # Model training\
  trainResults <- list()
  for(i in 1:length(CVindex)) {
    
    # Sample indices for repeated CV
    index <- list(CVindex[[i]])
    
    # Select samples from specific fold
    X_CV <- dataMatrix[,index[[1]]]
    factor_CV_cor <- ifelse(factor[index[[1]]] == "Yes",1,0)
    
    
    # Calculate correlations with factors
    correlations_CV <- apply(X_CV, 1, 
                             function(x){cor(x, factor_CV_cor,method = "spearman")})
    
    names(correlations_CV) <- rownames(X_CV)
    
    # Select top correlated features for each factor
    finalProbes <- names(tail(sort(abs(correlations_CV)),10000))
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = FALSE,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               index = index)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix[finalProbes,]),
                 y = factor(factor, levels = c("No", "Yes")),
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
save(outputList, file = "PerFactor/RandomForest/OutputList_factors_Cor_RF.RData")


################################################################################

# Model Performance

################################################################################
load("PerFactor/RandomForest/OutputList_factors_Cor_RF.RData")

PerformanceList <- list()
for (i in 1:length(outputList)){
  trainResults <- outputList[[i]]
  perf <- matrix(NA, nrow = 35, ncol = 25)
  for (j in 1:length(trainResults)){
    perf[,j] <- trainResults[[j]]$ROC
  }
  optPar <- which.max(rowMeans(perf))
  optPerf <- perf[optPar,]
  opt_mtry <- trainResults[[1]]$mtry[optPar]
  opt_splitrule <- trainResults[[1]]$splitrule[optPar]
  opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]
  
  PerformanceList[[i]] <- list(optPerf, opt_mtry, opt_splitrule, opt_min.node.size)
}
names(PerformanceList) <- names(outputList)
save(PerformanceList, file = "PerFactor/RandomForest/PerformanceList_factors_Cor_RF.RData")

################################################################################

# Obs vs pred in CV

################################################################################

# Load data
load("Data/X_test.RData")
load("Data/X_nonTest.RData")
load("PerFactor/Y_test_factors.RData")
load("PerFactor/Y_nonTest_factors.RData")
load("PerFactor/selectedProbes.RData")
selectedProbes <- selectedProbes[-13]

# Test if samples are in correct order
Y_train <- Y_nonTest[,-13] # ignore kidney disease
Y_test <- Y_test[,-13] # ignore kidney disease
X_train <- log2(X_nonTest/(1-X_nonTest))
X_test <- log2(X_test/(1-X_test))
all(colnames(X_train) == rownames(Y_train))
all(colnames(X_test) == rownames(Y_test))

load("PerFactor/RandomForest/PerformanceList_factors_Cor_RF.RData")

ObsPred_all <- list()
for (f in 1:length(PerformanceList)){
  
  # mtry
  mtry_CV <- PerformanceList[[f]][[2]]
  
  # splitrule
  splitrule_CV <- PerformanceList[[f]][[3]]
  
  # min.node.size
  min.node.size_CV <- PerformanceList[[f]][[4]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
  colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
  
  # ML method
  MLmethod = "ranger"
  
  # Performance metric
  performance_metric = "ROC"
  
  # factor
  factor <- Y_train[!is.na(Y_train[,f]),f]
  
  # dataMatrix
  dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
  
  # Create index
  set.seed(123)
  CVindex <- NULL
  for (r in 1:5){
    temp <- createFolds(1:length(factor),5, returnTrain = TRUE)
    names(temp) <- paste0(names(temp),".Rep",r)
    CVindex <- c(CVindex, temp)
  }
  rm(temp)
  
  # Model training
  ObsPred_CV <- NULL
  for(i in 1:length(CVindex)) {
    
    # Sample indices for repeated CV
    index <- list(CVindex[[i]])
    
    # Select samples from specific fold
    X_CV <- dataMatrix[,index[[1]]]
    factor_CV_cor <- ifelse(factor[index[[1]]] == "Yes",1,0)
    
    
    # Calculate correlations with factors
    correlations_CV <- apply(X_CV, 1, 
                             function(x){cor(x, factor_CV_cor,method = "spearman")})
    
    names(correlations_CV) <- rownames(X_CV)
    
    # Select top correlated features for each factor
    finalProbes <- names(tail(sort(abs(correlations_CV)),10000))
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(search = "grid", 
                               savePredictions = TRUE,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               index = index)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(dataMatrix[finalProbes,]),
                 y = factor(factor, levels = c("No", "Yes")),
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
save(ObsPred_all, file = "PerFactor/RandomForest/ObsPred_all_RF.RData")


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
save(finalOutput, file = "PerFactor/finalOutput_CV.RData")

factorNames <- c("Education", "Systolic BP", "BMI", "Total Chol.", "Physical Act.",
                 "Diet", "Smoking", "Alcohol Intake", "Depression", "Diabetes",
                 "HDL Chol.", "Heart Disease", "APOE", "Sex", "Age < 47", "Age > 53")

# Plot cross-validation performance
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

# Save plot
ggsave(plotCV, file = "PerFactor/Random Forest/CVperformance_RF_perFactor.png", width = 8, height = 6)




################################################################################

# Obs vs pred test set

################################################################################


# Load data
load("Data/X_test.RData")
load("Data/X_nonTest.RData")
load("PerFactor/Y_test_factors.RData")
load("PerFactor/Y_nonTest_factors.RData")
load("PerFactor/selectedProbes.RData")
load("PerFactor/PerformanceList_factors_Cor.RData")
load("PerFactor/finalOutput_CV.RData")
selectedProbes <- selectedProbes[-13]

# Test if samples are in correct order
Y_train <- Y_nonTest[,-13] # ignore kidney disease
Y_test <- Y_test[,-13] # ignore kidney disease
X_train <- log2(X_nonTest/(1-X_nonTest))
X_test <- log2(X_test/(1-X_test))
all(colnames(X_train) == rownames(Y_train))
all(colnames(X_test) == rownames(Y_test))


finalOutput_test <- list()
for (f in 1:length(PerformanceList)){
  
  # mtry
  mtry_CV <- PerformanceList[[f]][[2]]
  
  # splitrule
  splitrule_CV <- PerformanceList[[f]][[3]]
  
  # min.node.size
  min.node.size_CV <- PerformanceList[[f]][[4]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
  colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")
  
  # ML method
  MLmethod = "ranger"
  
  # Performance metric
  performance_metric = "ROC"
  
  # factor
  factor1 <- Y_train[!is.na(Y_train[,f]),f]
  factor_cor <- ifelse(factor1 == "Yes",1,0)
  
  # dataMatrix
  dataMatrix <- X_train[selectedProbes[[f]],!is.na(Y_train[,f])]
  
  # Calculate correlations with factors
  correlations_all <- apply(dataMatrix, 1, 
                            function(x){cor(x, factor_cor,method = "spearman")})
  
  names(correlations_all) <- rownames(dataMatrix)
  
  # Select top correlated features for each factor
  finalProbes <- names(tail(sort(abs(correlations_all)),10000))
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(method = "none", 
                             classProbs = TRUE)
  
  # Actual training
  set.seed(123)
  finalModel <- train(x = t(dataMatrix[finalProbes,]),
                      y = factor(factor1, levels = c("No", "Yes")),
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = TRUE)
  
  # Get test predictions
  pred <- predict(finalModel, t(X_test[finalProbes,]), type = "prob")
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
save(finalOutput_test, file = "PerFactor/RandomForest/finalOutput_test.RData")



# Performance in test set:

factorNames <- c("Education", "Systolic BP", "BMI", "Total Chol.", "Physical Act.",
                 "Diet", "Smoking", "Alcohol Intake", "Depression", "Diabetes",
                 "HDL Chol.", "Heart Disease", "APOE", "Sex", "Age < 47", "Age > 53")

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

# Save plot
ggsave(plotTest, file = "PerFactor/Random Forest/Testperformance_RF_perFactor.png", width = 8, height = 6)



