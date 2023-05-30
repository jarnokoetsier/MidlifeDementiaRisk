################################################################################

# Get best model for each factor

################################################################################

# All factors except for age, smoking, and sex
factors <- c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", "Physical", "HeartDisease",
             "Education", "Depression", "SysBP", "Diet", "SexMale", "Smoking")
factorNames <- c("BMI", "Type II Diabetes", "L-M Alchol Intake", "HDL Chol.",
                 "Total Chol.", "Physical Act.", "Heart Disease", "Education",
                 "Depression", "Systolic BP", "Dietary Intake", "Sex", "Smoking")

bestModel <- rep(NA, length(factors))
for (i in 1:length(factors)){
  
  if ((factors[i] != "SexMale")&(factors[i] != "Smoking")){
    # Literature-based feature selection + EN
    load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
    auc_CV_lit <- finalOutput[[factors[i]]]$AUC
    
    # Correlation-based feature selection + EN
    load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
    auc_CV_EN <- finalOutput[[factors[i]]]$AUC
    
    # Correlation-based feature selection + RF
    load("~/PerFactor/Random Forest/finalOutput_CV.RData")
    auc_CV_RF <- finalOutput[[factors[i]]]$AUC 
    
    bestModel[i] <- which.max(c(auc_CV_lit,auc_CV_EN, auc_CV_RF))
  }
  if (factors[i] == "SexMale"){
    bestModel[i] <- 2
  }
  if (factors[i] == "Smoking"){
    bestModel[i] <- 2
  }

}

bestModel
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



allModels <- list()
for (f in 1:length(factors)){
  
  if (bestModel[f] == 2){
    load("PerFactor/ElasticNet/PerformanceList_factors_Cor.RData")
    # Set grid for lambda
    lambdaCV <- PerformanceList[[factors[f]]][[3]]
    
    # Set grid for alpha
    alphaCV <- PerformanceList[[factors[f]]][[2]]
    
    # Combine into a single data frame
    parameterGrid <- expand.grid(alphaCV, lambdaCV)
    colnames(parameterGrid) <- c(".alpha", ".lambda")
    
    # ML method
    MLmethod = "glmnet"
  }

  # Performance metric
  performance_metric = "ROC"
  
  # factor
  factor1 <- Y_train[!is.na(Y_train[,factors[f]]),factors[f]]
  factor_cor <- ifelse(factor1 == "Yes",1,0)
  
  # dataMatrix
  dataMatrix <- X_train[selectedProbes[[factors[f]]],!is.na(Y_train[,factors[f]])]
  
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
  
  allModels[[f]] <- finalModel
}
factors <- c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", "Physical", "HeartDisease",
             "Education", "Depression", "SysBP", "Diet", "SexMale", "Smoking")
names(allModels) <- factors
save(allModels, file = "PerFactor/allModels.RData")
