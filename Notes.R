# Set working directory
setwd("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA")

# Load packages
library(caret)
library(glmnet)
library(foreach)
library(prospectr)

# Load data
load("E:/mydat1.RData")
load("E:/pheno.RData")

# Prepare data
all_Data <- as.matrix(t(mydat1))[!is.na(pheno$SmokingScore),]
all_Class <- pheno$SmokingScore[!is.na(pheno$SmokingScore)]

# Split data into training and test set
selectedSamples <- prospectr::kenStone(
  X = all_Data, 
  k = round(0.8*nrow(all_Data)),
  pc = 10,
  .center = TRUE,
  .scale = TRUE
)

# Make training and test data
trainingData <- all_Data[selectedSamples$model,]
trainingClass <- all_Class[selectedSamples$model]

testData <- all_Data[selectedSamples$test,]
testClass <- all_Class[selectedSamples$test]

# EN_regression: repeated cross-validation
EN_regression <- function(trainingData, trainingClass, nfold, nrep, alphaCV, lambdaCV){
  
   # Parameter grid: each combination of alpha and lambda
  parameterGrid <- expand.grid(alphaCV, lambdaCV)
  npar <- nrow(parameterGrid)
  
  # Empty matrix to collect model performance
  rmse <- matrix(NA, nrow = nfold*nrep, ncol = npar)
  mae <- matrix(NA, nrow = nfold*nrep, ncol = npar)
  
  count = 0
  # For each repeat.....
  for(r in 1:nrep){
    
    # Create folds
    folds <- createFolds(1:nrow(trainingData), k = nfold)
    
    for (i in 1:nfold){
      count = count + 1
      
      # Split training data in training and validation set:
      # Training
      X_train <- as.matrix(trainingData[-folds[[i]],])
      Y_train <- trainingClass[-folds[[i]]]
      
      # Validation
      X_val <- as.matrix(trainingData[folds[[i]],])
      Y_val <- trainingClass[folds[[i]]]
      
      for (p in 1:npar){
       
         # Fit model
        en_model_cv <- glmnet(x = X_train, 
                              y = Y_train, 
                              family = "gaussian",
                              alpha = parameterGrid[p,1], 
                              lambda = parameterGrid[p,1],
                              standardize = TRUE)
        
        # Get prediction
        pred <- predict(en_model_cv, X_val)
        
        # Use RMSE and MAE to evaluate model
        rmse[count,p] <- RMSE(pred = pred[,1], obs = Y_val)
        mae[count,p] <- MAE(pred = pred[,1], obs = Y_val)
      }
    }
  }
  
  # Save accuracies and parameters in a list object
  output <- list(rmse, mae, parameterGrid[,1], parameterGrid[,2])
  names(output) <- c("RMSE", "MAE", "Alpha", "Lambda")
  return(output)
}

nfold = 5
nrep = 1

# Alpha's and lambda's
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))
alphaCV <- seq(0,1,length.out = 10)

test <- EN_regression(trainingData, trainingClass, nfold, nrep, alphaCV, lambdaCV)








# Alpha's and lambda's
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))
alphaCV <- seq(0,1,length.out = 10)

# Parameter grid: each combination of alpha and lambda
parameterGrid <- expand.grid(alphaCV, lambdaCV)

# For each parameter combination....
output <- foreach(i =  1:nrow(parameterGrid), .packages = c("glmnet","caret"), .inorder = FALSE) %do% {
  
  # set seed
  set.seed(123)
  
  # Perform repeated cross-validation
  EN_regression(trainingData, trainingClass, nfold, nrep, parameterGrid[i,1], parameterGrid[i,2])
}

#fit = glmnet(as.matrix(mtcars[-1]), mtcars[,1], lambda=cv.glmnet(as.matrix(mtcars[-1]), mtcars[,1])$lambda.1se)


alphaCV <- seq(0,1,length.out = 11)
mse_1se <- rep(NA, length(alphaCV))
lambda <- rep(NA, length(alphaCV))
foldid = sample(rep(seq(5), length = nrow(trainingData)))
for (a in 1:length(alphaCV)){
  CVresults <- cv.glmnet(x = trainingData,
                      y = trainingClass,
                      foldid = foldid,
                      type.measure = "mse",
                      keep = FALSE,
                      alpha = alphaCV[a])
  mse_1se[a] <- CVresults$cvm[CVresults$index[2]]
  lambda[a] <- CVresults$lambda[CVresults$index[2]]
}

# Fit model
en_model <- glmnet(x = trainingData, 
                   y = trainingClass, 
                   alpha = alphaCV[which.min(mse_1se)], 
                   lambda = lambda[which.min(mse_1se)],
                   standardize = TRUE)


pred <- predict(en_model, trainingData)
RMSE(pred,trainingClass)
plot(pred, trainingClass)

caret::R2(pred,trainingClass)

pred <- predict(en_model, testData)
RMSE(pred,testClass)
plot(pred, testClass)

caret::R2(pred,testClass)

coefficients <- as.matrix(lambda$glmnet.fit$beta)



# annotation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
