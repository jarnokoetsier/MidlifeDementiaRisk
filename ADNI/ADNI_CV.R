library(caret)
library(glmnet)

# Clear workspace and console
rm(list = ls())
cat("\014") 


load("~/ADNI/Y_ADNI.RData")
load("~/ADNI/predictedScore.RData")

Y_train <- Y_ADNI[!is.na(Y_ADNI$conversion),]
X_train <- predictedScore[Y_train$Basename,-15]

#Create index
set.seed(123)
CVindex <- NULL
for (r in 1:5){
  temp <- createFolds(1:nrow(X_train),5, returnTrain = TRUE)
  names(temp) <- paste0(names(temp),".Rep",r)
  CVindex <- c(CVindex, temp)
}
rm(temp)

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE,
                           index = CVindex)


# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0.1,1,length.out = 10)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Use MSE as performance metric
performance_metric = "ROC"
MLmethod = "glmnet"

# Actual training
set.seed(123)
fit <- train(x = X_train,
             y = ifelse(Y_train$conversion == "converter", "Yes", "No"),
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)
save(fit, file = "ADNI_fit_EN_converter.RData")

set.seed(123)
fit <- train(x = X_train,
             y = ifelse(Y_train$TwoClass == 1, "Class1", "Class2"),
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = TRUE)
save(fit, file = "ADNI_fit_EN_TwoClass.RData")


# Get results
load("ADNI_fit_EN_converter.RData")
trainResults <- fit$results
optLambda <- fit$bestTune$lambda
optAlpha <- fit$bestTune$alpha
predictions <- fit$pred[(fit$pred$alpha == optAlpha) & (fit$pred$lambda == optLambda),]

boxplot(predictions$Yes ~ predictions$obs)