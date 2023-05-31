
# Load packages
library(tidyverse)
library(caret)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("~/EMIF/X_EMIF_imp.RData")
load("EMIF/metaData_EMIF.RData")
metaData_EMIF <- as.data.frame(metaData_EMIF)
rownames(metaData_EMIF) <- metaData_EMIF$X

# Remove individuals with unmatching sex
metaData_EMIF <- metaData_EMIF[metaData_EMIF$sex.match == 1,]

# Keep midlife samples only
metaData_EMIF <- metaData_EMIF[metaData_EMIF$Age <= 75,]

# Filter for CSF biomarkers only
rownames(metaData_EMIF ) <- metaData_EMIF $X
CSFbio <- metaData_EMIF [,c("Local_PTAU_Abnormal", "Local_TTAU_Abnormal", "Local_AB42_Abnormal")]
colnames(CSFbio) <- c("PTAU", "TTAU", "AB")
samples <- rownames(CSFbio)[(!is.na(CSFbio$PTAU)) & 
                              (!is.na(CSFbio$TTAU)) &
                              (!is.na(CSFbio$AB))]

metaData_CSF <- metaData_EMIF[samples,]
X_EMIF_fil <- X_EMIF_imp[rownames(X_EMIF_imp) %in% samples,]
metaData_CSF <- metaData_CSF[rownames(X_EMIF_fil),]
all(rownames(metaData_CSF) == rownames(X_EMIF_fil))

#*****************************************************************************#
# Split data
#*****************************************************************************#

# Male
nTrain0 <- round(0.7*nrow(X_EMIF_fil[metaData_CSF$Gender == 0,]))
selectedSamples0 <- prospectr::kenStone(
  X = X_EMIF_fil[metaData_CSF$Gender == 0,], 
  k = nTrain0,
  pc = 8,
  .center = TRUE,
  .scale = TRUE
)
temp0 <- rownames(X_EMIF_fil)[metaData_CSF$Gender == 0]
test0 <- temp0[selectedSamples0$test]
train0 <- temp0[selectedSamples0$model]

# Female
nTrain1 <- round(0.7*nrow(X_EMIF_fil[metaData_CSF$Gender == 1,]))
selectedSamples1 <- prospectr::kenStone(
  X = X_EMIF_fil[metaData_CSF$Gender == 1,], 
  k = nTrain1,
  pc = 8,
  .center = TRUE,
  .scale = TRUE
)
temp1 <- rownames(X_EMIF_fil)[metaData_CSF$Gender == 1]
test1 <- temp1[selectedSamples1$test]
train1 <- temp1[selectedSamples1$model]


# Combine into single train and test set
X_train <- X_EMIF_fil[c(train0, train1),]
Y_train <- metaData_CSF[c(train0, train1),]

X_test <- X_EMIF_fil[c(test0, test1),]
Y_test <- metaData_CSF[c(test0, test1),]

intersect(rownames(X_train), rownames(X_test))

sum(duplicated(X_test))
sum(duplicated(X_train))

table(Y_test$Diagnosis)
table(Y_train$Diagnosis)

# Save training and test set
save(X_train, file = "PredictCSF/X_train_PredictCSF.RData")
save(Y_train, file = "PredictCSF/Y_train_PredictCSF.RData")
save(X_test, file = "PredictCSF/X_test_PredictCSF.RData")
save(Y_test, file = "PredictCSF/Y_test_PredictCSF.RData")

#*****************************************************************************#
# correlation-based selection
#*****************************************************************************#

# Load training and test set
load("PredictCSF/X_train_PredictCSF.RData")
load("PredictCSF/Y_train_PredictCSF.RData")
load("PredictCSF/X_test_PredictCSF.RData")
load("PredictCSF/Y_test_PredictCSF.RData")

# Perform correlation-based selection
Y_train1 <- Y_train[,c("Local_PTAU_Abnormal", "Local_TTAU_Abnormal", "Local_AB42_Abnormal")]
correlations <- matrix(NA, nrow = ncol(X_train), ncol = ncol(Y_train1))
for (i in 1:ncol(Y_train1)) {
  correlations[,i] <- apply(X_train, 2, 
                            function(x){cor(x, 
                                            Y_train1[,i], 
                                            method = "spearman")})
}
rownames(correlations) <- colnames(X_train)
colnames(correlations) <- colnames(Y_train1)
save(correlations, file = "PredictCSF/correlations.RData")


# Selected probes
selectedProbes <- list()
for (i in 1:ncol(correlations)){
  selectedProbes[[i]] <- names(tail(sort(abs(correlations[,i])),70000))
}
names(selectedProbes) <- colnames(correlations)
save(selectedProbes, file = "PredictCSF/selectedProbes.RData")

################################################################################

# Model Training (ElasticNet)

################################################################################

# Load packages
library(tidyverse)
library(caret)
library(glmnet)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("PredictCSF/X_train_PredictCSF.RData")
load("PredictCSF/Y_train_PredictCSF.RData")
load("PredictCSF/X_test_PredictCSF.RData")
load("PredictCSF/Y_test_PredictCSF.RData")
load("PredictCSF/selectedProbes.RData")

Y_train1 <- Y_train[,c("Local_PTAU_Abnormal", "Local_TTAU_Abnormal", "Local_AB42_Abnormal")]

# Test if samples are in correct order
all(rownames(X_train) == rownames(Y_train1))
all(rownames(X_test) == rownames(Y_test))

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)


# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Machine learning method
MLmethod = "glmnet"
performance_metric = "ROC"

# Perform machine learning
CVList <- list()
for (f in 1:ncol(Y_train1)){
  
  # factor
  factor <- Y_train1[,f]
  
  # dataMatrix
  dataMatrix <- X_train[,selectedProbes[[f]]]
  
  # Model training\
  trainResults <- list()
  for(i in 1:length(CVindex)) {
    
    # Sample indices for repeated CV
    index <- list(CVindex[[i]])
    
    # Settings for repeated cross-valBasename
    fitControl <- trainControl(search = "grid", 
                               savePredictions = FALSE,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               index = index)
    
    
    # Select samples from specific fold
    X_CV <- dataMatrix[index[[1]],]
    factor_CV_cor <- factor[index[[1]]]
    
    # Calculate correlations with factors
    correlations_CV <- apply(X_CV, 2, 
                             function(x){cor(x, factor_CV_cor,method = "spearman")})
    
    names(correlations_CV) <- colnames(X_CV)
    
    # Select top correlated features for each factor
    finalProbes <- names(tail(sort(abs(correlations_CV)),10000))
    
    
    # Actual training
    set.seed(123)
    fit <- train(x = dataMatrix[,finalProbes],
                 y = ifelse(Y_train1[,f] == 1, "High", "Low"),
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)
    
    
    trainResults[[i]] <- fit$results
    
  }
  
  # Save output
  CVList[[f]] <- trainResults
  
  
}

names(CVList) <- colnames(Y_train1)
save(CVList, file = "PredictCSF/CVList.RData")

# Get optimal model performance and corresponding parameters in cross-validation
PerformanceList <- list()
for (i in 1:length(CVList)){
  trainResults <- CVList[[i]]
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
names(PerformanceList) <- names(CVList)
save(PerformanceList, file = "PredictCSF/PerformanceList.RData")


# Get best models
bestModels <- list()
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
  
  
  # dataMatrix
  dataMatrix <- X_train[,selectedProbes[[f]]]
  
  # Calculate correlations with factors
  correlations_all <- apply(dataMatrix, 2, 
                            function(x){cor(x, Y_train1[,f],method = "spearman")})
  
  names(correlations_all) <- colnames(dataMatrix)
  
  # Select top correlated features for each factor
  finalProbes <- names(tail(sort(abs(correlations_all)),10000))
  
  X_train_cor <- X_train[,finalProbes]
  X_test_cor <- X_test[,finalProbes]
  save(X_train_cor, X_test_cor, file = paste0("PredictCSF/X_all_",colnames(Y_train1)[f], ".RData"))
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(method = "none",
                             classProbs = TRUE)
  
  # Actual training
  set.seed(123)
  finalModel <- train(x = X_train_cor,
                      y = ifelse(Y_train1[,f] == 1, "High", "Low"),
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = TRUE)
  
  
  
  # Save
  bestModels[[f]] <- finalModel
  
}
names(bestModels) <- colnames(Y_train1)
save(bestModels, file = "PredictCSF/bestModels.RData")

# Evaluate models
factors <- c("Local_PTAU_Abnormal", "Local_TTAU_Abnormal", "Local_AB42_Abnormal")
f = 3
fit <- bestModels[[f]]
load(paste0("PredictCSF/X_all_", factors[f],".RData"))

test <- predict(fit, X_test_cor, type = "prob")
roc_list <- pROC::roc(Y_test[,factors[f]], test$High)
pROC::auc(roc_list)


################################################################################

# Model Training (ElasticNet, no correlation)

################################################################################

# Load packages
library(tidyverse)
library(caret)
library(glmnet)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("PredictCSF/X_train_PredictCSF.RData")
load("PredictCSF/Y_train_PredictCSF.RData")
load("PredictCSF/X_test_PredictCSF.RData")
load("PredictCSF/Y_test_PredictCSF.RData")
load("PredictCSF/selectedProbes.RData")

Y_train1 <- Y_train[,c("Local_PTAU_Abnormal", "Local_TTAU_Abnormal", "Local_AB42_Abnormal")]

# Test if samples are in correct order
all(rownames(X_train) == rownames(Y_train1))
all(rownames(X_test) == rownames(Y_test))

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)



# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Machine learning method
MLmethod = "glmnet"
performance_metric = "ROC"
fitControl <- trainControl(search = "grid", 
                           savePredictions = FALSE,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           index = CVindex)

# Model training
bestModels <- list()
for (i in 1:ncol(Y_train1)){
  set.seed(123)
  finalModel <- train(x = X_train,
                      y = ifelse(Y_train1[,f] == 1, "High", "Low"),
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = TRUE)
  
  
  bestModels[[f]] <- finalModel
}
names(bestModels) <- colnames(Y_train1)
save(bestModels, file = "PredictCSF/bestModels_ENnoCor.RData")

# Evalute best models
factors <- c("Local_PTAU_Abnormal", "Local_TTAU_Abnormal", "Local_AB42_Abnormal")
f = 3
fit <- bestModels[[f]]
load(paste0("PredictCSF/X_all_", factors[f],".RData"))

test <- predict(fit, X_test_cor, type = "prob")
roc_list <- pROC::roc(Y_test[,factors[f]], test$High)
pROC::auc(roc_list)


################################################################################

# Model Training (Random Forest)

################################################################################

# Load packages
library(tidyverse)
library(caret)
library(ranger)
library(e1071)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("PredictCSF/X_train_PredictCSF.RData")
load("PredictCSF/Y_train_PredictCSF.RData")
load("PredictCSF/X_test_PredictCSF.RData")
load("PredictCSF/Y_test_PredictCSF.RData")
load("PredictCSF/selectedProbes.RData")

Y_train1 <- Y_train[,c("Local_PTAU_Abnormal", "Local_TTAU_Abnormal", "Local_AB42_Abnormal")]

# Test if samples are in correct order
all(rownames(X_train) == rownames(Y_train1))
all(rownames(X_test) == rownames(Y_test))

# Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)


# Number of randomly selected predictors
mtry_CV <- c(100,500,1000,1200,1500,1700,2000,3000,4000)

# split rule
splitrule_CV <- "gini"

# minimal node size
min.node.size_CV = 3:30

# Combine into a single data frame
parameterGrid <- expand.grid(mtry_CV, splitrule_CV, min.node.size_CV)
colnames(parameterGrid) <- c(".mtry", ".splitrule", ".min.node.size")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "ranger"

CVList <- list()
for (f in 1:ncol(Y_train1)){
  
  # factor
  factor <- Y_train1[,f]
  
  # dataMatrix
  dataMatrix <- X_train[,selectedProbes[[f]]]
  
  # Model training\
  trainResults <- list()
  for(i in 1:length(CVindex)) {
    
    # Sample indices for repeated CV
    index <- list(CVindex[[i]])
    
    # Settings for repeated cross-valBasename
    fitControl <- trainControl(search = "grid", 
                               savePredictions = FALSE,
                               summaryFunction = twoClassSummary,
                               classProbs = TRUE,
                               index = index)
    
    
    # Select samples from specific fold
    X_CV <- dataMatrix[index[[1]],]
    factor_CV_cor <- factor[index[[1]]]
    
    # Calculate correlations with factors
    correlations_CV <- apply(X_CV, 2, 
                             function(x){cor(x, factor_CV_cor,method = "spearman")})
    
    names(correlations_CV) <- colnames(X_CV)
    
    # Select top correlated features for each factor
    finalProbes <- names(tail(sort(abs(correlations_CV)),10000))
    
    
    # Actual training
    set.seed(123)
    fit <- train(x = dataMatrix[,finalProbes],
                 y = ifelse(Y_train1[,f] == 1, "High", "Low"),
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)
    
    
    trainResults[[i]] <- fit$results
    
  }
  
  # Save output
  CVList[[f]] <- trainResults
  
  
}

names(CVList) <- colnames(Y_train1)
save(CVList, file = "PredictCSF/CVList_RF.RData")

# Get optimal performances and corresponding parameter values
PerformanceList <- list()
for (i in 1:length(CVList)){
  trainResults <- CVList[[i]]
  perf <- matrix(NA, nrow = 252, ncol = 25)
  for (j in 1:length(trainResults)){
    perf[,j] <- trainResults[[j]]$ROC
  }
  optPar <- which.max(rowMeans(perf))
  optPerf <- perf[optPar,]
  opt_mtry <- trainResults[[1]]$mtry[optPar]
  opt_min.node.size <- trainResults[[1]]$min.node.size[optPar]
  opt_splitrule <- trainResults[[1]]$splitrule[optPar]
  
  PerformanceList[[i]] <- list(optPerf, opt_mtry, opt_min.node.size, opt_splitrule)
}
names(PerformanceList) <- names(CVList)
save(PerformanceList, file = "PredictCSF/PerformanceList_RF.RData")

# Get best models
bestModels <- list()
for (f in 1:length(PerformanceList)){
  
  # Set grid for mtry
  opt_mtry <- PerformanceList[[f]][[2]]
  
  # Set grid for minimum node size
  opt_min.node.size <- PerformanceList[[f]][[3]]
  
  # Set grid for splitrule
  opt_splitrule <- PerformanceList[[f]][[4]]
  
  # Combine into a single data frame
  parameterGrid <- expand.grid(opt_mtry, opt_min.node.size, opt_splitrule)
  colnames(parameterGrid) <- c(".mtry", ".min.node.size", ".splitrule")
  
  # ML method
  MLmethod = "ranger"
  
  # Performance metric
  performance_metric = "ROC"
  
  # Load data
  load(paste0("PredictCSF/X_all_",colnames(Y_train1)[f], ".RData"))
  
  # Settings for repeated cross-validation
  fitControl <- trainControl(method = "none",
                             classProbs = TRUE)
  
  # Actual training
  set.seed(123)
  finalModel <- train(x = X_train_cor,
                      y = ifelse(Y_train1[,f] == 1, "High", "Low"),
                      metric= performance_metric,
                      method = MLmethod,
                      tuneGrid = parameterGrid,
                      trControl = fitControl,
                      maximize = TRUE)
  
  
  
  # Save
  bestModels[[f]] <- finalModel
  
}
names(bestModels) <- colnames(Y_train1)
save(bestModels, file = "PredictCSF/bestModels_RF.RData")

# Evaluate models
factors <- c("Local_PTAU_Abnormal", "Local_TTAU_Abnormal", "Local_AB42_Abnormal")
f = 1
fit <- bestModels[[f]]
load(paste0("PredictCSF/X_all_", factors[f],".RData"))

test <- predict(fit, X_test_cor, type = "prob")
roc_list <- pROC::roc(Y_test[,factors[f]], test$High)
pROC::auc(roc_list)



################################################################################

# Make plot

################################################################################

# Load pROC package
library(pROC)

# Get ROC curve for each model in test set
score <- c( "_ENnoCor","", "_RF")
scoreName <- c("ElasticNet","Correlation + ElasticNet", "Correlation + Random Forest")
factors <- c("Local_PTAU_Abnormal", "Local_TTAU_Abnormal", "Local_AB42_Abnormal")
factorName <- c('p-tau', "t-tau", "amyloid-beta")
plotDF <- as.data.frame(testDF)
ROCplot <- NULL
aucValue <- rep(NA, length(score)*length(factors))
ci_low <- rep(NA, length(score)*length(factors))
ci_high <- rep(NA, length(score)*length(factors))
for (i in 1:length(score)){
  
  # Load best model
  load(paste0("PredictCSF/bestModels",score[i],".RData"))
  load("PredictCSF/X_test_PredictCSF.RData")
  
  for (f in 1:length(factors)){
    fit <- bestModels[[f]]
    if (i == 1){
      testPred <- predict(fit, X_test, type = "prob")
    } else{
      load(paste0("PredictCSF/X_all_", factors[f],".RData"))
      testPred <- predict(fit, X_test_cor, type = "prob")
    }
    test <- pROC::roc(Y_test[,factors[f]],testPred$High, direction = "<")
    temp <- data.frame(Sensitivity = test$sensitivities,
                       Specificity = test$specificities,
                       Factor = rep(factorName[f],length(test$specificities)),
                       Model = rep(scoreName[i], length(test$specificities)))
    
    ROCplot <- rbind.data.frame(ROCplot, temp)
    aucValue[f + (i-1)*3] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
    ci_low[f + (i-1)*3] <- format(round(ci(test)[1],2),nsmall = 2)
    ci_high[f + (i-1)*3] <- format(round(ci(test)[3],2),nsmall = 2)
  }
  
  
}

plotAUC <- data.frame(AUC_ci = paste0("AUC: ",aucValue, " (", ci_low, "-", ci_high, ")"),
                      AUC = paste0("AUC: ",aucValue),
                      Factor = rep(factorName, 3),
                      Model = rep(scoreName, each = 3),
                      X = 0.9,
                      Y = rep(c(0.050,0.200,0.125), each = 3))



# Set colors
colorList <- list(c("#4292C6","#2171B5","#084594"),
                  c("#807DBA","#6A51A3","#4A1486"),
                  c("#EF3B2C","#CB181D", "#99000D")
)

# Make plot
n = 1 # Which CSF biomarker?
colors <- colorList[[n]]
p <- ggplot(ROCplot[ROCplot$Factor == factorName[n],]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Model), 
            size = 1.5, linetype = "solid") +
  geom_text(data = plotAUC[plotAUC$Factor == factorName[n],], aes(x = X, y = Y, label = AUC, color = Model),
            fontface = "bold") +
  scale_color_manual(values = colors) +
  ggtitle(paste0(factorName[n], " status")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = paste0("PredictCSF/ROC_",factorName[n],".png"), width = 7.5, height = 5)
