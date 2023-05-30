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
library(pROC)

# Load data
load("~/Data/X_nonTest.RData")
load("~/PerFactor/Y_nonTest_factors.RData")
load("~/Data/X_test.RData")
load("~/PerFactor/Y_test_factors.RData")


###############################################################################

# Smoking score

###############################################################################

# Smoking score function
smokingScore<-function(betas){
  
  load("Data/SmokingScoreRefData.rda")
  #contains Illig_data, Illig_data_up and Illig_data_down
  
  #subsetting own data 
  #SABRE_data_down
  select <- rownames(betas) %in% Illig_data_down$cpgs
  A_down <- subset(betas, select =="TRUE")
  
  #SABRE_data_up
  select <- rownames(betas) %in% Illig_data_up$cpgs
  A_up <- subset(betas, select =="TRUE")
  
  #sort SABRE data by Cpg name
  A_up <- A_up[order(rownames(A_up)),]
  A_down <- A_down[order(rownames(A_down)),]
  
  #match Illig data by by Cpg name
  Illig_data_up<-Illig_data_up[match(rownames(A_up), Illig_data_up$cpgs),]
  Illig_data_down<-Illig_data_down[match(rownames(A_down), Illig_data_down$cpgs),]
  
  #as some outliers have been removed and replaced with NAs need to handle missing values
  matrix_up_A<- matrix(nrow=nrow(Illig_data_up), ncol=ncol(A_up))
  for (i in 1:ncol(A_up)){
    matrix_up_A[,i]<- (A_up[,i])-(Illig_data_up$reference_never_median_beta_all)}
  colnames(matrix_up_A)<- colnames(A_up)
  rownames(matrix_up_A)<- Illig_data_up$cpgs
  
  #Calculate scores - ##UP##
  #calculate scores for each individual in the dataset
  scores_up_A<- as.numeric(rep(NA,ncol(A_up)))
  for (i in 1:ncol(A_up)){
    scores_up_A[i]<-sum(matrix_up_A[,i]*Illig_data_up$weights)}
  
  
  #Calculate diffs between SABRE beta values and the reference for each site - ##DOWN###
  matrix_down_A<- matrix(nrow=nrow(Illig_data_down), ncol=ncol(A_down))
  for (i in 1:ncol(A_down)){
    matrix_down_A[,i]<- (Illig_data_down$reference_never_median_beta_all)-(A_down[,i])}
  colnames(matrix_down_A)<- colnames(A_down)
  rownames(matrix_down_A)<- Illig_data_down$cpgs
  
  #Calculate scores - ##DOWN##
  #calculate scores for each individual in the dataset
  scores_down_A<- as.numeric(rep(NA,ncol(A_down)))
  for (i in 1:ncol(A_down)){
    scores_down_A[i]<-sum(matrix_down_A[,i]*Illig_data_down$weights)}
  
  ##combine scores
  scores_combined_A <- scores_up_A + scores_down_A
  
  return(scores_combined_A)
}

# Calculate smoking score
smoking_nonTest <- smokingScore(X_nonTest)
smoking_test <- smokingScore(X_test)

# Performance in training set:
roc_list <- roc(response = factor(Y_nonTest$Smoking,
                                  levels = c("No", "Yes")), 
                predictor = smoking_nonTest)

# AUC
AUC <- as.numeric(auc(roc_list))

AUCplot <- data.frame(Sensitivity = roc_list$sensitivities,
                      Specificity = roc_list$specificities)

# optimal threshold
Gmean <- sqrt(roc_list$sensitivities*roc_list$specificities)
threshold <- roc_list$thresholds[which.max(Gmean)]

# Observed versus predicted
ObsPred_CV <- data.frame(predictedClass = ifelse(predAge_nonTest[,i] < threshold, "Yes", "No"),
                         predictedAge = predAge_nonTest[,i],
                         obs = Y_nonTest$Smoking)

# Combine in list object
OutputSmoking_CV <- list(AUCplot = AUCplot, 
                         AUC = AUC, 
                         threshold = threshold, 
                         ObsPred_test = ObsPred_CV)

# Performance in test set
roc_list <- roc(response = factor(Y_test$Smoking,
                                  levels = c("No", "Yes")), 
                predictor = smoking_test)

# AUC
AUC <- as.numeric(auc(roc_list))

AUCplot <- data.frame(Sensitivity = roc_list$sensitivities,
                      Specificity = roc_list$specificities)

# optimal threshold
Gmean <- sqrt(roc_list$sensitivities*roc_list$specificities)
threshold <- roc_list$thresholds[which.max(Gmean)]

# Observed versus predicted
ObsPred_test <- data.frame(predictedClass = ifelse(predAge_test[,i] < threshold, "Yes", "No"),
                           predictedAge = predAge_test[,i],
                           obs = Y_test$Smoking)

# Combine in list object
OutputSmoking_test <- list(AUCplot = AUCplot, 
                           AUC = AUC, 
                           threshold = threshold, 
                           ObsPred_test = ObsPred_test)
  
  
# Save performances
save(OutputSmoking_CV, file = "PerFactor/OutputSmoking_CV.RData")
save(OutputSmoking_test, file = "PerFactor/OutputSmoking_test.RData")


#*****************************************************************************#
#   Make plots
#*****************************************************************************#

# AUC training
load("PerFactor/OutputSmoking_CV.RData")
auc_test_smoking <- OutputSmoking_test$AUC

load("~/PerFactor/ElasticNet/finalOutput_test.RData")
auc_test_smoking <- c(auc_test_smoking,finalOutput_test$Smoking$AUC)

load("~/PerFactor/Random Forest/finalOutput_test.RData")
auc_test_smoking <- c(auc_test_smoking,finalOutput_test$Smoking$AUC)

# AUC test
load("PerFactor/OutputSmoking_test.RData")
auc_CV_smoking <- OutputSmoking_CV$AUC

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
auc_CV_smoking <- c(auc_CV_smoking,finalOutput$Smoking$AUC)

load("~/PerFactor/Random Forest/finalOutput_CV.RData")
auc_CV_smoking <- c(auc_CV_smoking,finalOutput$Smoking$AUC)

# Prepare data for plotting
plotDF_smoking <- data.frame(AUC = c(auc_test_smoking,auc_CV_smoking),
                             Model = factor(rep(c("Smoking Score", "ElasticNet*", "Random Forest*"),2),
                                            levels = c("ElasticNet*", "Random Forest*","Smoking Score")),
                             Set = c(rep("Training",3), rep("Test",3))
)

# Set colors
color <- c("#1B9E77", "#D95F02", "#E7B10A")

# Make plot
p <- ggplot(plotDF_smoking) +
  geom_bar(aes(x = Model, y = AUC, fill = Model),
           stat = "identity", position = position_dodge(), color = "grey") +
  facet_grid(cols = vars(Set), scale = "free") +
  scale_fill_manual(values = color) +
  scale_alpha_manual(values = c(0.7,1)) +
  xlab(NULL) +
  ggtitle("Smoking Models") +
  guides(fill = "none") +
  coord_flip(ylim = c(0.7,1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "SmokingModels_AUC.png", width = 8, height = 3)

###############################################################################

# Age

###############################################################################

# Load packages
library(wateRmelon)
library(pROC)

# predicted age in training set
predAge_nonTest <- agep(X_nonTest, method='all')
predAge_nonTest <- predAge_nonTest[,c(1,3,5,7,9)]
save(predAge_nonTest, file = "PerFactor/predAge_nonTest.RData")

# predicted age in test set
predAge_test <- agep(X_test, method='all')
predAge_test <- predAge_test[,c(1,3,5,7,9)]
save(predAge_test, file = "PerFactor/predAge_test.RData")

# Load data
load("PerFactor/predAge_test.RData")
load("PerFactor/predAge_nonTest.RData")
modelName <- c("Horvath", "Hannum", "Pheno-Age", "Skin-Blood", "Lin")

# Get models' performances in training and test set
OutputAge47_CV <- list()
OutputAge53_CV <- list()
OutputAge47_test <- list()
OutputAge53_test <- list()

for (i in 1:ncol(predAge_nonTest)){
  
  # Age 47:
  
  # train
  roc_list <- roc(response = factor(Y_nonTest$Age47,
                                    levels = c("No", "Yes")), 
                  predictor = predAge_nonTest[,i])
  
  AUC <- as.numeric(auc(roc_list))
  
  AUCplot <- data.frame(Sensitivity = roc_list$sensitivities,
                     Specificity = roc_list$specificities,
                     Model = modelName[i])
  
  Gmean <- sqrt(roc_list$sensitivities*roc_list$specificities)
  threshold <- roc_list$thresholds[which.max(Gmean)]
  
  ObsPred_CV <- data.frame(predictedClass = ifelse(predAge_nonTest[,i] < threshold, "Yes", "No"),
                           predictedAge = predAge_nonTest[,i],
                           obs = Y_nonTest$Age47)
  
  OutputAge47_CV[[i]] <- list(AUCplot = AUCplot, 
                              AUC = AUC, 
                              threshold = threshold, 
                              ObsPred_test = ObsPred_CV)
  
  # test
  roc_list <- roc(response = factor(Y_test$Age47,
                                    levels = c("No", "Yes")), 
                  predictor = predAge_test[,i])
  
  AUC <- as.numeric(auc(roc_list))
  
  AUCplot <- data.frame(Sensitivity = roc_list$sensitivities,
                        Specificity = roc_list$specificities,
                        Model = modelName[i])
  
  ObsPred_test <- data.frame(predictedClass = ifelse(predAge_test[,i] < threshold, "Yes", "No"),
                           predictedAge = predAge_test[,i],
                           obs = Y_test$Age47)
  
  OutputAge47_test[[i]] <- list(AUCplot = AUCplot, 
                                AUC = AUC, 
                                threshold = threshold, 
                                ObsPred_test = ObsPred_test)
  
  
  
  # Age 53:
  
  # train
  roc_list <- roc(response = factor(Y_nonTest$Age53,
                                    levels = c("No", "Yes")), 
                  predictor = predAge_nonTest[,i])
  
  AUC <- as.numeric(auc(roc_list))
  
  AUCplot <- data.frame(Sensitivity = roc_list$sensitivities,
                        Specificity = roc_list$specificities,
                        Model = modelName[i])
  
  Gmean <- sqrt(roc_list$sensitivities*roc_list$specificities)
  threshold <- roc_list$thresholds[which.max(Gmean)]
  
  ObsPred_CV <- data.frame(predictedClass = ifelse(predAge_nonTest[,i] > threshold, "Yes", "No"),
                           predictedAge = predAge_nonTest[,i],
                           obs = Y_nonTest$Age53)
  
  OutputAge53_CV[[i]] <- list(AUCplot = AUCplot, 
                              AUC = AUC, 
                              threshold = threshold, 
                              ObsPred_test = ObsPred_CV)
  
  # test
  roc_list <- roc(response = factor(Y_test$Age53,
                                    levels = c("No", "Yes")), 
                  predictor = predAge_test[,i])
  
  AUC <- as.numeric(auc(roc_list))
  
  AUCplot <- data.frame(Sensitivity = roc_list$sensitivities,
                        Specificity = roc_list$specificities,
                        Model = modelName[i])
  
  Gmean <- sqrt(roc_list$sensitivities*roc_list$specificities)
  threshold <- roc_list$thresholds[which.max(Gmean)]
  
  ObsPred_test <- data.frame(predictedClass = ifelse(predAge_test[,i] > threshold, "Yes", "No"),
                           predictedAge = predAge_test[,i],
                           obs = Y_test$Age53)
  
  OutputAge53_test[[i]] <- list(AUCplot = AUCplot, 
                                AUC = AUC, 
                                threshold = threshold, 
                                ObsPred_test = ObsPred_test)
  
}

names(OutputAge47_CV) <- modelName
names(OutputAge53_CV) <- modelName
names(OutputAge47_test) <- modelName
names(OutputAge53_test) <- modelName

# Save performances
save(OutputAge47_CV, file = "PerFactor/OutputAge47_CV.RData")
save(OutputAge53_CV, file = "PerFactor/OutputAge53_CV.RData")
save(OutputAge47_test, file = "PerFactor/OutputAge47_test.RData")
save(OutputAge53_test, file = "PerFactor/OutputAge53_test.RData")


#*****************************************************************************#
#   Make plots
#*****************************************************************************#

# Prepare data for plotting
auc_CV_age47 <- rep(NA, length(modelName))
auc_CV_age53 <- rep(NA, length(modelName))
auc_test_age47 <- rep(NA, length(modelName))
auc_test_age53 <- rep(NA, length(modelName))

for (i in 1:length(OutputAge47_CV)){
  auc_CV_age47[i] <- OutputAge47_CV[[i]]$AUC
  auc_test_age47[i] <- OutputAge47_test[[i]]$AUC
  auc_CV_age53[i] <- OutputAge53_CV[[i]]$AUC
  auc_test_age53[i] <- OutputAge53_test[[i]]$AUC
}



load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
auc_CV_age47 <- c(auc_CV_age47,finalOutput$Age47$AUC)
auc_CV_age53 <- c(auc_CV_age53,finalOutput$Age53$AUC)


load("~/PerFactor/Random Forest/finalOutput_CV.RData")
auc_CV_age47 <- c(auc_CV_age47,finalOutput$Age47$AUC)
auc_CV_age53 <- c(auc_CV_age53,finalOutput$Age53$AUC)

load("~/PerFactor/ElasticNet/finalOutput_test.RData")
auc_test_age47 <- c(auc_test_age47,finalOutput$Age47$AUC)
auc_test_age53 <- c(auc_test_age53,finalOutput$Age53$AUC)

load("~/PerFactor/Random Forest/finalOutput_test.RData")
auc_test_age47 <- c(auc_test_age47,finalOutput$Age47$AUC)
auc_test_age53 <- c(auc_test_age53,finalOutput$Age53$AUC)


# AUC in cross-validation
plotAUC_CV <- data.frame(Model = c(modelName, c("ElasticNet*", "Random Forest*"),
                                modelName, c("ElasticNet*", "Random Forest*")),
                      Age = c(rep("Age < 47", length(auc_CV_age47)),
                              rep("Age > 53", length(auc_CV_age53))),
                      AUC = c(auc_CV_age47, auc_CV_age53)
)
plotAUC_CV$Set <- rep("Training", nrow(plotAUC_CV))

# AUC in cross-validation
plotAUC_test <- data.frame(Model = c(modelName, c("ElasticNet*", "Random Forest*"),
                                   modelName, c("ElasticNet*", "Random Forest*")),
                         Age = c(rep("Age < 47", length(auc_test_age47)),
                                 rep("Age > 53", length(auc_test_age53))),
                         AUC = c(auc_test_age47, auc_test_age53)
)
plotAUC_test$Set <- rep("Test", nrow(plotAUC_test))


# Combine AUC in test and cross-validation into single data frame
plotAUC <- rbind.data.frame(plotAUC_CV, plotAUC_test)
plotAUC$Model <- factor(plotAUC$Model,
                        levels = c("ElasticNet*", "Random Forest*",
                                   modelName))


# Set colors
color <- c(RColorBrewer::brewer.pal(n = 7, "Dark2"),
           RColorBrewer::brewer.pal(n = 7, "Set2"))

# Make plot
p <- ggplot() +
  geom_bar(data = plotAUC, aes(x = Model, y = AUC, alpha = Age, fill = Model),
           stat = "identity", position = position_dodge(), color = "grey") +
  facet_grid(cols = vars(Set), scale = "free") +
  scale_fill_manual(values = color) +
  scale_alpha_manual(values = c(0.7,1)) +
  xlab(NULL) +
  ggtitle("Age Models") +
  guides(fill = "none") +
  coord_flip(ylim = c(0.85,1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "AgeModels_AUC.png", width = 8, height = 6)



# Set colors
color <- c(RColorBrewer::brewer.pal(n = 5, "Dark2"),
           RColorBrewer::brewer.pal(n = 5, "Set2"))

# Prepare data
plotDF <- rbind.data.frame(plotDF_age47, plotDF_age53)
plotDF$Age <- c(rep("Age < 47",nrow(plotDF_age47)),
                rep("Age > 53",nrow(plotDF_age53)))
plotDF$Group <- paste0(plotDF$Age, ", ", plotDF$Model)

# Make ROC curves
p <- ggplot() +
  geom_path(data = plotDF, aes(y = Sensitivity, x = 1- Specificity,
                               color = Group), 
            linewidth = 1, linetype = "solid") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  scale_color_manual(values = color) +
  ggtitle("Age") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))