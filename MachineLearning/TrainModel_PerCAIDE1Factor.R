# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(glmnet)
library(spls)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)

# Load machine learning functions
source("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA/MachineLearning/FUN_MachineLearning.R")

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load data
load("cellType.RData")
load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

Y_CAIDE1_factors$Sex <- ifelse(Y_CAIDE1_factors$Sex == 1, "Male", "Female")
Y_CAIDE1_factors$Edu_c <- ifelse(Y_CAIDE1_factors$Edu_c == 0, "Educated", "Uneducated")
Y_CAIDE1_factors$PHYSICAL_c <- ifelse(Y_CAIDE1_factors$PHYSICAL_c == 0, "Active", "Inactive")

#=============================================================================#
# FILL IN
#=============================================================================#

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "var"

# Load data
files <- list.files(paste0("X_", FeatureSelection))
for (f in files){
  load(paste0("X_", FeatureSelection, "/", f))
}

# Prepare data
X_train = log2(X_CAIDE1_varMCor/(1-X_CAIDE1_varMCor))
#X_train = t(X_CAIDE1_PC)
Y_train = Y_CAIDE1_factors[,8:14]

# Test if samples are in correct order
all(colnames(X_train) == Y_CAIDE1_factors$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5

#=============================================================================#



###############################################################################

# ElasticNet

###############################################################################


#*****************************************************************************#
# Model training
#*****************************************************************************#

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
for (variables in 1:ncol(Y_train)){
  
  
  #****************************************************************************#
  # Classification
  #****************************************************************************#
  if (class(Y_train[,variables]) == "character"){
    
    # Performance metric
    performance_metric = "Kappa"
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "repeatedcv", 
                               number = nfold, 
                               repeats = nrep, 
                               search = "grid", 
                               savePredictions = FALSE)
    
    # Register cores for parallel computing
    nCores <- 3
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(X_train),
                 y = Y_train[,variables],
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)
    
    # Stop clusters
    stopCluster(cl)
    
    # Get results
    trainResults <- fit$results
    
    # Get optimal lambda and alpha
    optAlpha <- fit$bestTune$alpha
    optLambda <- fit$bestTune$lambda
    
    # Get coefficients, prediction, and performance during the repeated CV
    coefs <- matrix(NA, ncol(t(X_train))+1, nfold*nrep)
    perf <- rep(NA, nfold*nrep)
    folds <- fit$control$index
    pred_CV_class <- NULL
    pred_CV_response <- NULL
    obs_CV <- NULL
    fold_CV <- NULL
    count = 0
    for (r in 1:nrep){
      for (f in 1:nfold){
        count = count + 1
        en_model_cv <- glmnet(x = t(X_train)[folds[[count]],], 
                              y = Y_train[folds[[count]] ,variables], 
                              family = "binomial",
                              alpha = optAlpha, 
                              lambda = optLambda,
                              standardize = TRUE)
        
        # Get coefficients
        coefs[,count] <- as.matrix(coef(en_model_cv))
        
        # Get prediction
        pred_class <- predict(en_model_cv, t(X_train)[-folds[[count]],], type = "class")[,1]
        pred_response <- predict(en_model_cv, t(X_train)[-folds[[count]],], type = "response")[,1]
        
        # Get performance
        perf[count] <- sum(pred_class == Y_train[-folds[[count]],variables])/length(pred_class)
        
        pred_CV_class <- c(pred_CV_class,pred_class)
        pred_CV_response <- c(pred_CV_response,pred_response)
        obs_CV <- c(obs_CV, Y_train[-folds[[count]],variables])
        fold_CV <- c(fold_CV, rep(count,length(pred_class)))
      }
    }
    
    # Save observed and predicted in a dataframe
    ObsPred_CV <- data.frame(Predicted_class = pred_CV_class,
                             Predicted_response = pred_CV_response,
                             Observed = obs_CV,
                             Fold = fold_CV)
    
    # Get final model
    finalModel <- glmnet(x = t(X_train), 
                         y = Y_train[,variables], 
                         family = "binomial",
                         alpha = optAlpha, 
                         lambda = optLambda,
                         standardize = TRUE)
    
    # Save output
    outputList[[variables]] <- list(trainResults, 
                                    optLambda, 
                                    optAlpha, 
                                    perf, 
                                    ObsPred_CV, 
                                    coefs, 
                                    finalModel)

  } 
  
  #****************************************************************************#
  # Regression
  #****************************************************************************#
  if (class(Y_train[,variables]) != "character"){
    # Performance metric
    performance_metric = "RMSE"
    
    # Settings for repeated cross-validation
    fitControl <- trainControl(method = "repeatedcv", 
                               number = nfold, 
                               repeats = nrep, 
                               search = "grid", 
                               savePredictions = FALSE,
                               summaryFunction = regressionSummary)
    
    # Register cores for parallel computing
    nCores <- 3
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    # Actual training
    set.seed(123)
    fit <- train(x = t(X_train),
                 y = Y_train[,variables],
                 metric= performance_metric,
                 method = MLmethod,
                 tuneGrid = parameterGrid,
                 trControl = fitControl,
                 maximize = TRUE)
    
    # Stop clusters
    stopCluster(cl)
    
    # Get results
    trainResults <- fit$results
    
    # Get optimal lambda and alpha
    optAlpha <- fit$bestTune$alpha
    optLambda <- fit$bestTune$lambda
    
    # Get coefficients, prediction, and performance during the repeated CV
    coefs <- matrix(NA, ncol(t(X_train))+1, nfold*nrep)
    perf <- rep(NA, nfold*nrep)
    folds <- fit$control$index
    pred_CV <- NULL
    obs_CV <- NULL
    fold_CV <- NULL
    count = 0
    for (r in 1:nrep){
      for (f in 1:nfold){
        count = count + 1
        en_model_cv <- glmnet(x = t(X_train)[folds[[count]],], 
                              y = Y_train[folds[[count]],variables], 
                              family = "gaussian",
                              alpha = optAlpha, 
                              lambda = optLambda,
                              standardize = TRUE)
        
        # Get coefficients
        coefs[,count] <- as.matrix(coef(en_model_cv))
        
        # Get prediction
        pred <- predict(en_model_cv, t(X_train)[-folds[[count]],])[,1]
        
        # Get performance
        perf[count] <- RMSE(pred = pred, obs = Y_train[-folds[[count]],variables])
        
        pred_CV <- c(pred_CV,pred)
        obs_CV <- c(obs_CV, Y_train[-folds[[count]],variables])
        fold_CV <- c(fold_CV, rep(count,length(pred)))
      }
    }
    
    # Save observed and predicted in a dataframe
    ObsPred_CV <- data.frame(Predicted = pred_CV,
                             Observed = obs_CV,
                             Fold = fold_CV)
    
    # Get final model
    finalModel <- glmnet(x = t(X_train), 
                         y = Y_train[,variables], 
                         family = "gaussian",
                         alpha = optAlpha, 
                         lambda = optLambda,
                         standardize = TRUE)
    
    # Save output
    outputList[[variables]] <- list(trainResults, 
                                    optLambda, 
                                    optAlpha, 
                                    perf, 
                                    ObsPred_CV, 
                                    coefs, 
                                    finalModel)
    
  }


}

names(outputList) <- colnames(Y_train)
save(outputList, file = paste0("OutputList_CAIDE1factors_", FeatureSelection, ".RData"))


###############################################################################

# Evaluate

###############################################################################

methods = c("S", "var", "varM", "varCor", "varMCor", "PC")
methodsName <- c("S-score", "Variance (\u03b2)", "Variance (M)", 
                 "Variance (\u03b2, Cor)", "Variance (M, Cor)", "PCA")
plot_continuous <- NULL
plot_discrete <- NULL
for (m in 1:length(methods)){
  load(paste0("PerFactor/OutputList_CAIDE1factors_", methods[m], ".RData"))
  
  # Continuous variables
  continuous <- c("Age", "MeanSysBP", "BMI", "Chol_unloged")
  Rsquared <- NULL
  for (i in continuous){
    ObsPred_CV <- outputList[[i]][[5]]
    temp <- R2(pred = ObsPred_CV$Predicted, obs = ObsPred_CV$Observed)
    Rsquared <- c(Rsquared,temp)
  }
  
  temp1 <- data.frame(Variable =  c("Age", "Blood Pressure", "BMI", "Cholesterol"),
                     R2 = Rsquared,
                     Method = rep(methodsName[m], length(continuous)))
  
  plot_continuous <- rbind.data.frame(plot_continuous, temp1)
  
  # Discrete variables
  discrete <- c("Sex", "Edu_c", "PHYSICAL_c")
  kappa <- NULL
  for (i in discrete){
    ObsPred_CV <- outputList[[i]][[5]]
    temp <- confusionMatrix(factor(ObsPred_CV$Predicted_class), factor(ObsPred_CV$Observed))
    kappa <- c(kappa,temp$overall[2])
  }
  temp2 <- data.frame(Variable =  c("Sex", "Education", "Physical"),
                     Kappa = kappa,
                     Method = rep(methodsName[m], length(discrete)))
  
  plot_discrete <- rbind.data.frame(plot_discrete, temp2)
}


plot_continuous$Method <- factor(plot_continuous$Method, 
                                 levels = methodsName)


p_continuous <- ggplot(plot_continuous) +
  geom_bar(aes(x = Method, y = R2, fill = Method), 
           stat = "identity", color = "black") +
  facet_grid(rows = vars(Variable)) +
  xlab("Feature Selection Method") +
  ylab(expression(R^2)) +
  ylim(c(0,1)) +
  ggtitle("Continuous Variables") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "grey", linewidth = 0, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())

plot_discrete$Method <- factor(plot_discrete$Method, 
                                 levels = methodsName)


p_discrete <- ggplot(plot_discrete) +
  geom_bar(aes(x = Method, y = Kappa, fill = Method), 
           stat = "identity", color = "black") +
  facet_grid(rows = vars(Variable)) +
  xlab("Feature Selection Method") +
  ylab("Cohen's \u03BA") +
  ggtitle("Discrete Variables") +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.background = element_rect(fill= "grey", linewidth = 0, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))


library(patchwork)

p <- p_continuous / p_discrete +
  plot_layout(heights = c(5,4))

ggsave(p, file = "perCAIDE1factor_performance.png", width = 9, height = 9)



###############################################################################

# Predict CAIDE1

###############################################################################

methods = c("S", "var", "varM", "varCor", "varMCor", "PC")
methodsName <- c("S-score", "Variance (\u03b2)", "Variance (M)", 
                 "Variance (\u03b2, Cor)", "Variance (M, Cor)", "PCA")
predList <- list()
obsList <- list()
output <- NULL
output_obs <- NULL
for (m in 1:length(methods)){
  load(paste0("PerFactor/OutputList_CAIDE1factors_", methods[m], ".RData"))
  
  # Continuous variables
  continuous <- c("Age", "MeanSysBP", "BMI", "Chol_unloged")
  for (i in continuous){
    ObsPred_CV <- outputList[[i]][[5]]
    temp <- ObsPred_CV$Predicted
    temp_obs <- ObsPred_CV$Observed
    if (!is.null(output)){
      output <- cbind.data.frame(output,temp)
      output_obs <- cbind.data.frame(output_obs,temp_obs)
    } else{
      output <- temp
      output_obs <- temp_obs
    }
  }
  # Discrete variables
  discrete <- c("Sex", "Edu_c", "PHYSICAL_c")
  for (i in discrete){
    ObsPred_CV <- outputList[[i]][[5]]
    
    temp <- ObsPred_CV$Predicted_class
    temp_obs <- ObsPred_CV$Observed
    
    output <- cbind.data.frame(output,temp)
    output_obs <- cbind.data.frame(output_obs,temp_obs)
  }
  colnames(output) <- c(continuous, discrete)
  output$Fold <- ObsPred_CV$Fold
  
  colnames(output_obs) <- c(continuous, discrete)
  output_obs$Fold <- ObsPred_CV$Fold
  
  predList[[m]] <- output
  obsList[[m]] <- output_obs
}
names(predList) <- methods
names(obsList) <- methods

# Calculate CAIDE score
calculateCAIDE1 <- function(x){
  score <- as.data.frame(matrix(NA, nrow = nrow(x), ncol = 7))
  colnames(score) <- c("Age", "Sex", "Edu", "BP", "BMI", "Chol", "Phys")
  
  # Age
  score$Age[x$Age < 47] <- 0
  score$Age[x$Age >= 47 & x$Age <= 53] <- 3
  score$Age[x$Age > 53] <- 4
  
  # Sex
  score$Sex <- ifelse(x$Sex == "Male",1,0)
  
  # Education
  score$Edu <- ifelse(x$Edu_c == "Uneducated",2,0)
  
  # Systolic blood pressure
  score$BP <- ifelse(x$MeanSysBP <= 140,0,2)
  
  # BMI
  score$BMI <-ifelse(x$BMI <= 30,0,2)
  
  # Serum total cholesterol
  score$Chol <- ifelse(x$Chol_unloged <= 6.5,0,2)
  
  # Physical activity
  score$Phys <- ifelse(x$PHYSICAL_c == "Active",0,1)
  
  return(rowSums(score))
}

plotDF <- NULL
for (k in 1:length(predList)){
  for (fold in 1:25){
    pred <- calculateCAIDE1(predList[[k]][predList[[k]]$Fold == fold,])
    obs <- calculateCAIDE1(obsList[[k]][obsList[[k]]$Fold == fold,])
    
    method <- names(predList)[k]
    rmse <- caret::RMSE(pred, obs)
    
    temp <- data.frame(Method = method,
                       Fold = fold,
                       RMSE = rmse)
    plotDF <- rbind.data.frame(plotDF, temp)
  }
}


# Score and feature selection method
Score = "CAIDE1"
methods = c("S", "var", "varM", "varCor", "varMCor", "PC")

# Retrieve performance in CV for the different feature selection methods
Performances <- list()
for (i in 1:length(methods)){
  
  # load output
  load(paste0("CV_", Score, "_", methods[i],".RData"))
  
  # Put performance into list
  Performances[[i]] <- perf
}
names(Performances) <- methods

# Combine into single data frame
performanceDF <- data.frame(
  Method = c(rep("S-score", 25),
              rep("Variance (\u03b2)", 25),
              rep("Variance (M)", 25),
              rep("Variance (\u03b2, Cor)",25),
              rep("Variance (M, Cor)", 25),
              rep("PCA", 25)),
  Fold = rep(1:25,length(methods)),
  RMSE = unlist(Performances)
)

plotDF <- rbind.data.frame(plotDF, performanceDF)
plotDF$Approach <- c(rep("Per Factor", 150), rep("Per Score", 150))

ggplot(plotDF) +
  geom_boxplot(aes(x = Method, y = RMSE, fill = Approach))




