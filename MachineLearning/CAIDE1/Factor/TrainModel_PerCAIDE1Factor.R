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

Y <- Y_CAIDE1[,c(1,2,5,15:20)]
Y$Sex_c <- ifelse(Y$Sex == 1, "Male", "Female")
Y$Edu_c <- ifelse(Y$Edu_c == 0, "Educated", "Uneducated")
Y$Syst_c <- ifelse(Y$Syst_c == 0, "Low", "High")
Y$BMI_c <- ifelse(Y$BMI_c == 0, "Low", "High")
Y$Chol_c <- ifelse(Y$Chol_c == 0, "Low", "High")
Y$PHYSICAL_c <- ifelse(Y$PHYSICAL_c == 0, "Active", "Inactive")

#=============================================================================#
# FILL IN
#=============================================================================#

# Score and feature selection method
Score = "CAIDE1"
FeatureSelection = "Cor"

# Load data
files <- list.files(paste0("X/X_", FeatureSelection))
for (f in files){
  load(paste0("X/X_", FeatureSelection, "/", f))
}

# Prepare data
X_train = log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))
Y_train = Y[,3:9]

# Test if samples are in correct order
all(colnames(X_train) == Y$Basename)

# Set number of folds and repeats
nfold = 5
nrep = 5

load("CVindex.RData")
#=============================================================================#

###############################################################################

# KS + Correlation

###############################################################################


###############################################################################

# Correlation

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
    performance_metric = "ROC"
    
    # Y_train
    Y_train <- Y[,variables]
    
    # Register cores for parallel computing
    nCores <- 3
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    trainResults <- foreach::foreach(i = 1:length(CVindex), .packages = c("caret", "glmnet")) %dopar% {
      
      # Select samples from specific fold
      index <- list(CVindex[[i]])
      X_CV <- X_train[,index[[1]]]
      
      # Calculate correlations (using X_CV)
      factors <- Y_CAIDE1[index[[1]],14:20]
      correlations_CV <- matrix(NA, nrow = nrow(X_CV), ncol = ncol(factors))
      for (f in 1:ncol(factors)) {
        correlations_CV[,f] <- apply(X_CV, 1, 
                                     function(x){cor(x, 
                                                     factors[,f], 
                                                     method = "spearman")})
      }
      rownames(correlations_CV) <- rownames(X_CV)
      colnames(correlations_CV) <- colnames(factors)
      
      # Select top correlated features for each factor
      probes <- list()
      for (p in 1:7){
        probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
      }
      
      # get exactly 10,000 probes
      n = 1
      finalProbes <- unique(unlist(probes))
      while (length(finalProbes) > 10000){
        probes[[n]] <- probes[[n]][-1]
        finalProbes <- unique(unlist(probes))
        
        if (n < 7){
          n = n + 1
        } else {
          n = 1
        }
      }
      
      # Settings for repeated cross-validation
      fitControl <- trainControl(search = "grid", 
                                 savePredictions = FALSE,
                                 summaryFunction = twoClassSummary,
                                 classProbs = TRUE,
                                 index = index)
      
      # Actual training
      set.seed(123)
      fit <- train(x = t(X_train[finalProbes,]),
                   y = Y_train,
                   metric= performance_metric,
                   method = MLmethod,
                   tuneGrid = parameterGrid,
                   trControl = fitControl,
                   maximize = TRUE)
      
      
      return(fit$results)
      
    }
    
    # Stop clusters
    stopCluster(cl)
    
    
    # Save output
    outputList[[variables]] <- trainResults

  } 
  
  #****************************************************************************#
  # Regression
  #****************************************************************************#
  if (class(Y_train[,variables]) != "character"){
    # Performance metric
    performance_metric = "RMSE"
    
    # Y_train
    Y_train <- Y[,variables]
    
    # Register cores for parallel computing
    nCores <- 3
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    
    trainResults <- foreach::foreach(i = 1:length(CVindex), .packages = c("caret", "glmnet")) %dopar% {
      
      # Select samples from specific fold
      index <- list(CVindex[[i]])
      X_CV <- X_train[,index[[1]]]
      
      # Calculate correlations (using X_CV)
      factors <- Y_CAIDE1[index[[1]],14:20]
      correlations_CV <- matrix(NA, nrow = nrow(X_CV), ncol = ncol(factors))
      for (f in 1:ncol(factors)) {
        correlations_CV[,f] <- apply(X_CV, 1, 
                                     function(x){cor(x, 
                                                     factors[,f], 
                                                     method = "spearman")})
      }
      rownames(correlations_CV) <- rownames(X_CV)
      colnames(correlations_CV) <- colnames(factors)
      
      # Select top correlated features for each factor
      probes <- list()
      for (p in 1:7){
        probes[[p]] <- names(tail(sort(abs(correlations_CV[,p])),1700))
      }
      
      # get exactly 10,000 probes
      n = 1
      finalProbes <- unique(unlist(probes))
      while (length(finalProbes) > 10000){
        probes[[n]] <- probes[[n]][-1]
        finalProbes <- unique(unlist(probes))
        
        if (n < 7){
          n = n + 1
        } else {
          n = 1
        }
      }
      
      # Settings for repeated cross-validation
      fitControl <- trainControl(search = "grid", 
                                 savePredictions = FALSE,
                                 summaryFunction = regressionSummary,
                                 classProbs = TRUE,
                                 index = index)
      
      # Actual training
      set.seed(123)
      fit <- train(x = t(X_train[finalProbes,]),
                   y = Y_train,
                   metric= performance_metric,
                   method = MLmethod,
                   tuneGrid = parameterGrid,
                   trControl = fitControl,
                   maximize = TRUE)
      
      
      return(fit$results)
      
    }
    
    # Stop clusters
    stopCluster(cl)
    
    
    # Save output
    outputList[[variables]] <- trainResults
  }
}

names(outputList) <- colnames(Y[,3:9])
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

for (m in 1:length(methods)){
  load(paste0("PerFactor/OutputList_CAIDE1factors_", methods[m], ".RData"))
  output <- NULL
  output_obs <- NULL
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
names(predList) <- methodsName
names(obsList) <- methodsName

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
plotDF$Method <- factor(plotDF$Method,
                        levels = methodsName)

colors <- RColorBrewer::brewer.pal(3, "Reds")
p <- ggplot(plotDF) +
  geom_boxplot(aes(x = Method, y = RMSE, fill = Approach), alpha = 0.3) +
  geom_point(aes(x = Method, y = RMSE, color = Approach), 
             position=position_jitterdodge(jitter.width = 0.3)) +
  scale_fill_manual(values = colors[2:3]) +
  scale_color_manual(values = colors[2:3]) +
  xlab("Feature Selection Method") +
  ggtitle("CAIDE1", subtitle = "Performance in the cross-validation") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, file= "CAIDE1_RMSE_boxplot_perFactor.png", width = 8, height = 5)



# R2 instead of RMSE

plotDF <- NULL
for (k in 1:length(predList)){
  for (fold in 1:25){
    pred <- calculateCAIDE1(predList[[k]][predList[[k]]$Fold == fold,])
    obs <- calculateCAIDE1(obsList[[k]][obsList[[k]]$Fold == fold,])
    
    
    method <- names(predList)[k]
    rmse <- caret::R2(pred, obs)
    
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
  
  perf <- NULL
  for (j in 1:25){
    pred <- ObsPred_CV$Predicted[ObsPred_CV$Fold == j]
    obs <- ObsPred_CV$Observed[ObsPred_CV$Fold == j]
    rsq <- caret::R2(pred, obs)
    perf <- c(perf, rsq)
  }
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
plotDF$Method <- factor(plotDF$Method,
                        levels = methodsName)

colors <- RColorBrewer::brewer.pal(3, "Reds")
p <- ggplot(plotDF) +
  geom_boxplot(aes(x = Method, y = RMSE, fill = Approach), alpha = 0.3) +
  geom_point(aes(x = Method, y = RMSE, color = Approach), 
             position=position_jitterdodge(jitter.width = 0.3)) +
  scale_fill_manual(values = colors[2:3]) +
  scale_color_manual(values = colors[2:3]) +
  xlab("Feature Selection Method") +
  ylab(expression(R^2)) +
  ggtitle("CAIDE1", subtitle = "Performance in the cross-validation") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     face = "italic",
                                     size = 10))

ggsave(p, file= "CAIDE1_RMSE_boxplot_perFactor_R2.png", width = 8, height = 5)

