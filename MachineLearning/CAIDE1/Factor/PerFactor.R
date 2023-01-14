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

# Make training and test data
load("CAIDE.Rdata")
load("CAIDE2.Rdata")
load("EPILIBRA.Rdata")
load("metaData_ageFil.RData")
load("methSet_allNorm_fil.RData")
load("cellType.RData")
load("TestTrain.RData")
load("CVindex_CAIDE1.RData")


#################################################################################

# Prepare methylation data (X)

################################################################################

# Remove X and Y chromosomal genes
load("probe_annotation.RData")
probe_ann_fil <- probe_annotation[(probe_annotation$Chr != "chrX") &
                                    (probe_annotation$Chr != "chrY"), ]
methSet_allNorm_fil <- methSet_allNorm_fil[rownames(methSet_allNorm_fil) %in% probe_ann_fil$ID,]

# Split data into training and test set
X_all <- as.matrix(t(methSet_allNorm_fil))
X_test <- t(X_all[TestTrain$Test,])
X_train <- t(X_all[setdiff(dat$Basename, TestTrain$Test),])

save(X_test, file = "X_test.RData")
save(X_train, file = "X_train.RData")

#################################################################################

# Prepare phenotype data (Y)

################################################################################

# High Education
Education <- ifelse(dat$None.of.the.above == 1,"Yes","No")

# High Systolic Blood Pressure
SysBP <- ifelse(dat$MeanSysBP <= 140,"No","Yes")

# High BMI
BMI <-ifelse(dat$BMI <= 30,"No","Yes")

# High total Cholesterol
TotalChol <- ifelse(as.numeric(dat$Chol_unloged) <= 6.5,"No","Yes")

# Low physical activity
Physical <- ifelse(dat$Exercise.increased.pulse.more.than.2halfhrsawk == 1,"No","Yes")

# Healthy Diet
Diet <- ifelse(as.numeric(dat$Fruit) > 3 | as.numeric(dat$Vegtables) > 3, "Yes", "No")

# Smoking
Smoking <- ifelse(dat$Do.you.currently.smoke==1,"Yes","No")

# Low alcohol intake
Alcohol <- rep(NA,nrow(dat))
Alcohol[dat$How.often.do.you.drink.alcohol==0 ,"LtoMAlcohol"] <- "Yes"
Alcohol[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==1]<- "Yes"
Alcohol[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==1]<-"Yes"
Alcohol[dat$alcoholic.drinks.per.day==3 &    dat$How.often.do.you.drink.alcohol==1]<-"Yes"
Alcohol[dat$alcoholic.drinks.per.day==4 &    dat$How.often.do.you.drink.alcohol==1]<-"Yes"
Alcohol[dat$alcoholic.drinks.per.day==5 &    dat$How.often.do.you.drink.alcohol==1]<-"Yes"

Alcohol[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==2]<-"Yes"
Alcohol[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==2]<-"Yes"

Alcohol[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==3]<-"Yes"
Alcohol[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==3]<-"Yes"

# High alcohol intake
Alcohol[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==1]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==3  &  dat$How.often.do.you.drink.alcohol==2]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==2]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==2]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==3]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==3]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==3]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==2 &   dat$How.often.do.you.drink.alcohol==4]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==4]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==4]<-"No"
Alcohol[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==4]<-"No"
Alcohol[dat$alcoholic.drinks.per.day=="." &    dat$How.often.do.you.drink.alcohol==4]<-"No"

# Depression
Depression <- ifelse(dat$Depression==1,"Yes","No")

# Diabetes
Diabetes <-ifelse(dat$T2.Diabetes==1,"Yes","No")

# High HDL
HDL <- ifelse(dat$HDL_unloged > 2.2, "Yes","No")

# Heart Disease
HeartDisease <- ifelse(dat$Heart.Disease==1,"Yes","No")

# Kidney Disease
KidneyDisease <- ifelse(dat$Kidney..Disease==1,"Yes","No")


Y_all <- data.frame(Education,
                    SysBP,
                    BMI,
                    TotalChol,
                    Physical,
                    Diet,
                    Smoking,
                    Alcohol,
                    Depression,
                    Diabetes,
                    HDL,
                    HeartDisease,
                    KidneyDisease
)

rownames(Y_all) <- dat$Basename

Y_test <- Y_all[colnames(X_test),] 
Y_train <- Y_all[colnames(X_train),] 

save(Y_test, file = "Y_test.RData")
save(Y_train, file = "Y_train.RData")

#################################################################################

# Prepare correlation data

################################################################################

# Load data
load("X_test.RData")
load("X_train.RData")
load("Y_test.RData")
load("Y_train.RData")

all(colnames(X_train) == rownames(Y_train))
all(colnames(X_test) ==  rownames(Y_test))

#*****************************************************************************#
# correlation-based selection
#*****************************************************************************#

correlations <- matrix(NA, nrow = nrow(X_train), ncol = ncol(Y_train))
for (i in 1:ncol(Y_train)) {
  factor <- Y_train[!is.na(Y_train[,i]),i]
  factor <- ifelse(factor == "Yes",1,0)
  dataMatrix <- X_train[,!is.na(Y_train[,i])]
  
  correlations[,i] <- apply(dataMatrix, 1, 
                            function(x){cor(x, 
                                            factor, 
                                            method = "spearman")})
}
rownames(correlations) <- rownames(X_test)
colnames(correlations) <- colnames(Y_train)
save(correlations, file = "correlations_factors.RData")


# Selected probes
selectedProbes <- list()
for (i in 1:ncol(correlations)){
  selectedProbes[[i]] <- names(tail(sort(abs(correlations[,i])),50000))
}
names(selectedProbes) <- colnames(correlations)
save(selectedProbes, file = "selectedProbes.RData")

################################################################################

# Model Training

################################################################################

# Load data
load("X_test.RData")
load("X_train.RData")
load("Y_test.RData")
load("Y_train.RData")
load("selectedProbes.RData")

# Test if samples are in correct order
X_train <- log2(X_train/(1-X_train))
X_test <- log2(X_test/(1-X_test))
all(colnames(X_train) == rownames(Y_train))
all(colnames(X_test) == rownames(Y_test))

Y_train <- Y_train[,1:12]
Y_test <- Y_test[,1:12]

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
                 y = factor,
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

names(outputList) <- colnames(Y_train)[1:12]
save(outputList, file = "OutputList_factors_Cor.RData")


################################################################################

# Model Performance

################################################################################
load("OutputList_factors_Cor.RData")
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
save(PerformanceList, file = "PerformanceList_factors_Cor.RData")

################################################################################

# Obs vs pred in CV

################################################################################


# Load data
load("X_train.RData")
load("Y_train.RData")
load("selectedProbes.RData")
load("PerformanceList_factors_Cor.RData")

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
                 y = factor,
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
save(ObsPred_all, file = "ObsPred_all.RData")



# construct ROC
for (f in 1:length(ObsPred_all)){
  ObsPred_CV <- ObsPred_all[[f]]
  roc_list_CV <- roc(response = factor(ObsPred_CV$Observed, levels = c("Yes", "No")), 
                     predictor = ObsPred_CV$Predicted)
  
  # Get sensitivies and specificities
  plotDF_lowRisk <- data.frame(Sensitivity = roc_list_CV$sensitivities,
                               Specificity = roc_list_CV$specificities)
  
  # Get AUC
  auc_lowRisk <- auc(roc_list_CV)
  
  # Get predictions
  Gmean_CV <- sqrt(roc_list_CV$sensitivities*roc_list_CV$specificities)
  threshold <- roc_list_CV$thresholds[which.max(Gmean_CV)]
  ObsPred_CV <- ObsPred_CV
  ObsPred_CV$PredictedClass <- ifelse(ObsPred_CV$Predicted > threshold, "Yes", "No")
}

# Performance on test set
pred <- predict(finalModel, t(testData), type = "response")
ObsPred_test <- data.frame(Observed = ifelse(Y_test$CAIDE < 4, 2, 1),
                           Predicted = as.numeric(pred))

roc_list_test <- roc(response = ObsPred_test$Observed, 
                     predictor = ObsPred_test$Predicted)


# Get sensitivies and specificities
plotDF_lowRisk_test <- data.frame(Sensitivity = roc_list_test$sensitivities,
                                  Specificity = roc_list_test$specificities)

# Get AUC
auc_lowRisk_test <- auc(roc_list_test)

ObsPred_test_lowRisk <- ObsPred_test
ObsPred_test_lowRisk$PredictedClass <- ifelse(ObsPred_test_lowRisk$Predicted > threshold, "Low", "High/Intermediate")
ObsPred_test_lowRisk$ObservedScore <- Y_test$CAIDE

