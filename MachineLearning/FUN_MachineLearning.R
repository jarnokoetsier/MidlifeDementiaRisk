

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

