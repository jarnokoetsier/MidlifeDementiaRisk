nfold = 5
nrep = 1
library(doParallel)

test <- function (pred, obs){
  isNA <- is.na(pred)
  pred <- pred[!isNA]
  obs <- obs[!isNA]
  if (!is.factor(obs) && is.numeric(obs)) {
    if (length(obs) + length(pred) == 0) {
      out <- rep(NA, 3)
    }
    else {
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
      mse <- mean((pred - obs)^2)
      mae <- mean(abs(pred - obs))
      out <- c(mse, sqrt(mse), resamplCor^2, mae)
    }
    names(out) <- c("MSE","RMSE", "Rsquared", "MAE")
  }
  else {
    if (length(obs) + length(pred) == 0) {
      out <- rep(NA, 2)
    }
    else {
      pred <- factor(pred, levels = levels(obs))
      requireNamespaceQuietStop("e1071")
      out <- unlist(e1071::classAgreement(table(obs, pred)))[c("diag", 
                                                               "kappa")]
    }
    names(out) <- c("Accuracy", "Kappa")
  }
  if (any(is.nan(out))) 
    out[is.nan(out)] <- NA
  out
}

mse <- function (data, lev = NULL, model = NULL){
  if (is.character(data$obs)) 
    data$obs <- factor(data$obs, levels = lev)
  test(data[, "pred"], data[, "obs"])
}

fitControl <- trainControl(method = "repeatedcv", 
                           number = nfold, 
                           repeats = nrep, 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = mse)

# Actual training
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))
alphaCV <- seq(0,1,length.out = 11)
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")
s_metric = "MSE"
MLmethod = "glmnet"

nCores <- detectCores() - 1
cl <- makeCluster(nCores)
registerDoParallel(cl)
fit <- train(x = trainingData,
             y = trainingClass,
             metric= s_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

stopCluster(cl)


trainResults <- fit$results



pred_train <- predict(fit, trainingData)
RMSE(pred = pred_train, obs = trainingClass)
plot(trainingClass, pred_train)

pred_test <- predict(fit, testData)
RMSE(pred = pred_test, obs = testClass)
plot(testClass, pred_test)