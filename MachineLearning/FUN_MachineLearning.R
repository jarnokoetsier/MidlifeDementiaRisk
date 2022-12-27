

#*****************************************************************************#
# regressionSummary
#*****************************************************************************#

regressionSummary <- function (data, lev = NULL, model = NULL){
  
  # Remove missing values from prediction
  isNA <- is.na(data$pred)
  pred <- data$pred[!isNA]
  obs <- data$obs[!isNA]

  # Calculate correlation coefficient (r)
  if (length(unique(pred)) < 2 || length(unique(obs)) < 
      2) {
    resamplCor <- NA
  }
  else {
    resamplCor <- try(cor(pred, obs, use = "pairwise.complete.obs"), 
                      silent = TRUE)
    if (inherits(resamplCor, "try-error")) 
      resamplCor <- NA
  }
  
  # Calculate Mean Squared Error (MSE)
  mse <- mean((pred - obs)^2)
  
  # Calculate Mean Absolute error (MAE)
  mae <- mean(abs(pred - obs))
  
  # Save MSE, RMSE, R2, and MAE as output
  out <- c(mse, sqrt(mse), resamplCor^2, mae)
  
  names(out) <- c("MSE","RMSE", "Rsquared", "MAE")
  return(out)
}


#*****************************************************************************#
# classificationSummary
#*****************************************************************************#

classificationSummary <- function (data, lev = NULL, model = NULL){
  if (length(lev) > 2) {
    stop(paste("Your outcome has", length(lev), "levels. The twoClassSummary() function isn't appropriate."))
  }
  caret::requireNamespaceQuietStop("pROC")
  if (!all(levels(data[, "pred"]) == lev)) {
    stop("levels of observed and predicted data do not match")
  }
  rocObject <- try(pROC::roc(data$obs, data[, lev[1]], direction = ">", 
                             quiet = TRUE), silent = TRUE)
  rocAUC <- if (inherits(rocObject, "try-error")) 
    NA
  else rocObject$auc
  
  sens <- sensitivity(data[, "pred"], data[, "obs"], lev[1])
  spec <-  specificity(data[, "pred"], data[, "obs"], lev[2])
  out <- c(rocAUC, 
           sens, 
           spec,
           sqrt(sens*spec))
  names(out) <- c("ROC", "Sens", "Spec", "GMean")
  out
}

#*****************************************************************************#
# selectionKS
#*****************************************************************************#

# Select features using a Kennard-Stone like algorithm:

# dataMatrix: M-values
# nFeatures: Number of features to select
# seedProbe: which probe should the algorithm start with

selectionKS <- function(dataMatrix, nFeatures = 5000, seedProbe = NULL) {
  
  # Make matrix to save correlations
  cor_matrix <- matrix(NA, nrow(dataMatrix), nFeatures-1)
  rownames(cor_matrix) <- rownames(dataMatrix)
  
  # Make vector to save feature set
  newProbe <- rep(NA, nFeatures)
  newProbe[1] <- seedProbe
  
  # Start selecting features
  for (j in 1:(nFeatures - 1)){
    
    # Calculate correlations between seed probe and all other probes
    cor_matrix[,j] <- apply(dataMatrix,1,function(x){cor(x,dataMatrix[newProbe[j],], method = "spearman")})
    
    # Add most uncorrelated probe to feature set
    newProbe[j+1] <- rownames(cor_matrix)[which.min(abs(cor_matrix[,j]))]
    
    # Remove all highly correlated probes: keep probes with cor < 0.5
    cor_matrix <- cor_matrix[abs(cor_matrix[,j]) < 0.5,]
    dataMatrix <- dataMatrix[rownames(cor_matrix),]
  }
  return(newProbe)
}

