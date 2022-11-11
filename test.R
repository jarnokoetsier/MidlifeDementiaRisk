
# Set data and output directory
dataDir <- "E:/Thesis/MLData"
outputDir <- dataDir

# Load packages
library(caret)
library(glmnet)
library(foreach)
library(prospectr)
library(doParallel)
library(tidyverse)
library(ggrepel)

#*****************************************************************************#
# Data Preparation
#*****************************************************************************#

# Load data
load(paste0(dataDir, "/mydat1.RData"))
load(paste0(dataDir, "/pheno.RData"))

# Prepare data
all_X <- as.matrix(t(mydat1))[!is.na(pheno$SmokingScore),]
all_Y <- pheno$SmokingScore[!is.na(pheno$SmokingScore)]

# Split data into training and test set
selectedSamples <- prospectr::kenStone(
  X = all_X, 
  k = round(0.8*nrow(all_X)),
  pc = 10,
  .center = TRUE,
  .scale = TRUE
)

# Make training and test data
X_train <- t(all_X[selectedSamples$model,])
Y_train <- all_Y[selectedSamples$model]

X_test <- t(all_X[selectedSamples$test,])
Y_test <- all_Y[selectedSamples$test]

# Save data
save(X_train, Y_train, X_test, Y_test, file = paste0(dataDir, "/splitData.RData"))


#*****************************************************************************#
# Check test and training data in PCA plot
#*****************************************************************************#

# Scale training data
X_train_scaled <- (X_train - rowMeans(X_train))/(apply(X_train,1,sd))

# Make PCA model
pcaList <-  prcomp(t(X_train_scaled),        
                   retx = TRUE,
                   center = FALSE,
                   scale = FALSE,
                   rank. = 10)

# Get the PCA scores of the training data
scores_train <- as.data.frame(pcaList$x)
scores_train$ID <- rownames(scores_train)
scores_train <- inner_join(scores_train, pheno, by = c("ID" = "X"))

# Calculate the explained variance of the PCs
explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)

# Scale test data (using the standard deviation and mean of the training data)
X_test_scaled <- t((X_test - rowMeans(X_train))/(apply(X_train,1,sd)))

# Calculate the scores of the test data
scores_test <- as.data.frame(as.matrix(X_test_scaled) %*% as.matrix(pcaList$rotation))
scores_test$ID <- rownames(scores_test)
scores_test <- inner_join(scores_test, pheno, by = c("ID" = "X"))

# Combine the scores of test and training data in a single data frame
scores_all <- rbind.data.frame(scores_train, scores_test)
scores_all$Train <- c(rep("Training", nrow(scores_train)), rep("Test", nrow(scores_test)))

# Plot the scores of the training data and the project scores of the test data 
PCA_TrainTest <- ggplot()+
  geom_point(data = scores_all, aes(x = PC1, y = PC2, shape = Train, color = Train), size = 2, alpha = 0.9) +
  scale_color_brewer(palette = "Set1") +
  xlab(paste0("PC1 (", explVar[1], "%)")) +
  ylab(paste0("PC2 (", explVar[2], "%)")) +
  labs(title = "PCA Score Plot", 
       caption = "NOTE: The PCA model is constructed using the training data only. The test data is projected.") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"))


#*****************************************************************************#
# Feature selection
#*****************************************************************************#




#*****************************************************************************#
# Repeated Cross-Validation
#*****************************************************************************#
gc()
rm(list = ls())
setwd("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA")
source("FUN_MachineLearning.R")

# Set data and output directory
dataDir <- "E:/Thesis/MLData"
outputDir <- dataDir

# Load data
load(paste0(dataDir, "/splitData.RData"))

# Set number of folds and repeats
nfold = 5
nrep = 10

# Settings for repeated cross-validation
fitControl <- trainControl(method = "repeatedcv", 
                           number = nfold, 
                           repeats = nrep, 
                           search = "grid", 
                           savePredictions = TRUE,
                           summaryFunction = regressionSummary)

# Set grid for lambda
lambdaCV <- exp(seq(log(0.01),log(2.5),length.out = 100))

# Set grid for alpha
alphaCV <- seq(0,1,length.out = 11)

# Combine into a single data frame
parameterGrid <- expand.grid(alphaCV, lambdaCV)
colnames(parameterGrid) <- c(".alpha", ".lambda")

# Use MSE as performance metric
performance_metric = "MSE"
MLmethod = "glmnet"

# Register cores for parallel computing
nCores <- detectCores() - 1
cl <- makeCluster(nCores)
registerDoParallel(cl)

# Actual training
fit <- train(x = X_train,
             y = Y_train,
             metric= performance_metric,
             method = MLmethod,
             tuneGrid = parameterGrid,
             trControl = fitControl,
             maximize = FALSE)

# Stop clusters
stopCluster(cl)

# Save model
save(fit, file = paste0(dataDir,"/fit.RData"))

trainResults <- fit$results

# Evaluate CV performance
ggplot() +
  geom_tile(data = trainResults, 
            aes(x = log(lambda), y = alpha, fill = RMSE)) +
  geom_point(data = trainResults[which.min(trainResults$MSE),], 
             aes(x = log(lambda), y = alpha), 
             color = 'red', size = 3) +
  geom_text_repel(data = trainResults[which.min(trainResults$MSE),], 
                   aes(x = log(lambda), y = alpha, label = paste0("Alpha: ", round(alpha,2), "\n",
                                                                  "Lambda: ", round(lambda,2), "\n",
                                                                  "Mean RMSE: ", round(RMSE,2), "\n",
                                                                  "SD RMSE: ", round(RMSESD,2))), 
                   color = "white", alpha = 0.9) +
  xlab("log Lambda") +
  ylab("Alpha") +
  ggtitle("RMSE in Repeated Cross-Validation") +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(legend.title = element_text(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"))
    

optAlpha <- trainResults$alpha[which.min(trainResults$MSE)]
optLambda <- trainResults$lambda[which.min(trainResults$MSE)]

test <- as.data.frame(as.matrix(coef(fit$finalModel)))
test <- colSums(abs(test))

coefs <- matrix(NA, ncol(X_train)+1, nfold*nrep)
count = 0
for (r in 1:nrep){
  folds <- createFolds(1:nrow(X_train), k = nfold)
  for (f in 1:nfold){
    count = count + 1
    en_model_cv <- glmnet(x = X_train[-folds[[f]],], 
                          y = Y_train[-folds[[f]]], 
                          family = "gaussian",
                          alpha = optAlpha, 
                          lambda = optLambda,
                          standardize = TRUE)
    coefs[,count] <- as.matrix(coef(en_model_cv))
  }
}

coefs1 <- coefs[-1,]
coefs1 <- coefs[rowSums(coefs) != 0,]
coefPlot <- gather(as.data.frame(coefs1))
coefPlot$cpg <- rep(1:nrow(coefs1), ncol(coefs1))

ggplot(coefPlot, aes(x = cpg, y = key, fill = value)) +
  geom_tile()+
  scale_fill_viridis_c()

test <- varImp(fit$finalModel)



minMSE <- min(trainResults$MSE)
minMSE_SE <- trainResults$MSESD[which.min(trainResults$MSE)]/sqrt(nrep*nfold)

lambda_1se <- max(trainResults$lambda[trainResults$MSE < minMSE + minMSE_SE])
lambda_min <- trainResults$lambda[which.min(trainResults$MSE)]


pred_train <- predict(fit, trainingData)
RMSE(pred = pred_train, obs = trainingClass)
plot(trainingClass, pred_train)

pred_test <- predict(fit, testData)
RMSE(pred = pred_test, obs = testClass)
plot(testClass, pred_test)