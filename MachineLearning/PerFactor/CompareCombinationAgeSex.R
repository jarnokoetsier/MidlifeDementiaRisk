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

setwd("E:/Thesis/EXTEND/Methylation")

################################################################################

# Get best model for each factor

################################################################################

# All factors except for age, smoking, and sex
factors <- c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", "Physical", "HeartDisease",
             "Education", "Depression", "SysBP", "Diet")
factorNames <- c("BMI", "Type II Diabetes", "L-M Alchol Intake", "HDL Chol.",
                 "Total Chol.", "Physical Inact.", "Heart Disease", "Education",
                 "Depression", "Systolic BP", "Dietary Intake")

bestModel_test <- list()
bestModel_CV <- list()
for (i in 1:length(factors)){
  # Literature-based feature selection + EN
  load("PerFactor/Literature/finalOutput_lit_CV.RData")
  auc_CV_lit <- finalOutput[[factors[i]]]$AUC
  
  # Correlation-based feature selection + EN
  load("PerFactor/ElasticNet/finalOutput_CV.RData")
  auc_CV_EN <- finalOutput[[factors[i]]]$AUC
  
  # Correlation-based feature selection + RF
  load("PerFactor/Random Forest/finalOutput_CV.RData")
  auc_CV_RF <- finalOutput[[factors[i]]]$AUC 
  
  best <- which.max(c(auc_CV_lit,auc_CV_EN, auc_CV_RF))
  if (best == 1){
    load("~/PerFactor/Literature/finalOutput_lit_test.RData")
    load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
  }
  if (best == 2){
    load("PerFactor/ElasticNet/finalOutput_test.RData")
    load("PerFactor/ElasticNet/finalOutput_CV.RData")
  }
  if (best == 3){
    load("PerFactor/Random Forest/finalOutput_test.RData")
    load("PerFactor/Random Forest/finalOutput_CV.RData")
  }
  bestModel_test[[i]] <- finalOutput_test[[factors[i]]]
  bestModel_CV[[i]] <- finalOutput[[factors[i]]]
}

# Sex
load("PerFactor/ElasticNet/finalOutput_test.RData")
bestModel_test[[12]] <- finalOutput_test[["SexMale"]]
load("PerFactor/ElasticNet/finalOutput_CV.RData")
bestModel_CV[[12]] <- finalOutput[["SexMale"]]

# Age < 47
load("PerFactor/OutputAge47_test.RData")
bestModel_test[[13]] <- OutputAge47_test$`Skin-Blood`
load("PerFactor/OutputAge47_CV.RData")
bestModel_CV[[13]] <- OutputAge47_CV$`Skin-Blood`

# Age > 53
load("PerFactor/OutputAge53_test.RData")
bestModel_test[[14]] <- OutputAge53_test$`Skin-Blood`
load("PerFactor/OutputAge53_CV.RData")
bestModel_CV[[14]] <- OutputAge53_CV$`Skin-Blood`

# Smoking
load("PerFactor/ElasticNet/finalOutput_test.RData")
bestModel_test[[15]] <- finalOutput_test$Smoking

load("PerFactor/ElasticNet/finalOutput_CV.RData")
bestModel_CV[[15]] <- finalOutput$Smoking


names(bestModel_test) <- c(factors, "SexMale","Age47", "Age53", "Smoking")
names(bestModel_CV) <- c(factors, "SexMale","Age47", "Age53", "Smoking")

################################################################################

# CAIDE1

################################################################################

#******************************************************************************#
# Continuous
#******************************************************************************#
load("Y/Y_test.RData")
perf_all <- rep(NA, 5)
perf_agesex <- rep(NA, 5)
perf_no_agesex <- rep(NA, 5)

# Age
Age47 <- bestModel_test$Age47$ObsPred_test$predictedClass
Age53 <- bestModel_test$Age53$ObsPred_test$predictedClass
Age <- rep(3,length(Age47))
Age[Age47 == "Yes"] <- 0
Age[Age53 == "Yes"] <- 4

# Sex
Sex <- ifelse(bestModel_test$SexMale$ObsPred_test$predClass == "Yes",1,0)

# BMI
threshold <- bestModel_test$BMI$threshold
power <- log(0.5)/log(threshold)
BMI <- (1-bestModel_test$BMI$ObsPred_test$pred^power)*2

# Education
threshold <- bestModel_test$Education$threshold
power <- log(0.5)/log(threshold)
Education <- (1-bestModel_test$Education$ObsPred_test$pred^power)*2

# Physical Activity
threshold <- bestModel_test$Physical$threshold
power <- log(0.5)/log(threshold)
Physical <- (1-bestModel_test$Physical$ObsPred_test$pred^power)*1

# Systolic blood pressure
threshold <- bestModel_test$SysBP$threshold
power <- log(0.5)/log(threshold)
SysBP <- (1-bestModel_test$SysBP$ObsPred_test$pred^power)*2

# Total Cholesterol
threshold <- bestModel_test$TotalChol$threshold
power <- log(0.5)/log(threshold)
TotalChol <- (1-bestModel_test$TotalChol$ObsPred_test$pred^power)*2

# Predicted score with all factors
predCAIDE1 <- Age + BMI + Education + Physical + Sex + SysBP + TotalChol
perf_all[1] <- R2(pred = predCAIDE1, obs = Y_test$CAIDE)

predCAIDE1 <- Age + Sex + mean(BMI) + mean(Education) + mean(Physical) + mean(SysBP) + mean(TotalChol)
perf_agesex[1] <- R2(pred = predCAIDE1, obs = Y_test$CAIDE)

predCAIDE1 <- mean(Age) + BMI + Education + Physical + mean(Sex) + SysBP + TotalChol
perf_no_agesex[1] <- R2(pred = predCAIDE1, obs = Y_test$CAIDE)


#******************************************************************************#
# Categorical
#******************************************************************************#

# Age
Age47 <- bestModel_test$Age47$ObsPred_test$predictedClass
Age53 <- bestModel_test$Age53$ObsPred_test$predictedClass
Age <- rep(3,length(Age47))
Age[Age47 == "Yes"] <- 0
Age[Age53 == "Yes"] <- 4

# Sex
Sex <- ifelse(bestModel_test$SexMale$ObsPred_test$predClass == "Yes",1,0)

# BMI
BMI <- ifelse(bestModel_test$BMI$ObsPred_test$predClass == "Yes",2,0)

# Education
Education <- ifelse(bestModel_test$Education$ObsPred_test$predClass == "Yes",2,0)

# Physical activity
Physical <- ifelse(bestModel_test$Physical$ObsPred_test$predClass == "Yes",1,0)

# Systolic Blood pressure
SysBP <- ifelse(bestModel_test$SysBP$ObsPred_test$predClass == "Yes",2,0)

# Total Cholesterol
TotalChol <- ifelse(bestModel_test$TotalChol$ObsPred_test$predClass == "Yes",2,0)


# Predicted score with all factors
predCAIDE1 <- Age + BMI + Education + Physical + Sex + SysBP + TotalChol
perf_all[2] <- R2(pred = predCAIDE1, obs = Y_test$CAIDE)

predCAIDE1 <- Age + Sex + mean(BMI) + mean(Education) + mean(Physical) + mean(SysBP) + mean(TotalChol)
perf_agesex[2] <- R2(pred = predCAIDE1, obs = Y_test$CAIDE)

predCAIDE1 <- mean(Age) + BMI + Education + Physical + mean(Sex) + SysBP + TotalChol
perf_no_agesex[2] <- R2(pred = predCAIDE1, obs = Y_test$CAIDE)

#******************************************************************************#
# Data-driven
#******************************************************************************#

methods <- c("EN", "sPLS", "RF")
load("PerFactor/CombineFactors/predictedScore_factors_EXTEND.RData")
for (i in 1:length(methods)){
  # Load model
  load(paste0("PerFactor/CombineFactors/Fit_CombineFactors_", "CAIDE1", "_", methods[i], ".RData"))
  
  
  # make prediction
  X_test <- predictedScore_factors[Y_test$Basename,]
  
  pred <- predict(fit, X_test)
  perf_all[i + 2] <- R2(pred = pred, obs = Y_test$CAIDE)
  

  X_test1 <- X_test
  for (j in 1:ncol(X_test)){
    X_test1[,j] <- rep(mean(X_test[,j]),nrow(X_test))
  }
  X_test1$Age <- X_test$Age 
  X_test1$SexMale <- X_test$SexMale
  pred <- predict(fit, X_test1)
  perf_agesex[i + 2] <- R2(pred = pred, obs = Y_test$CAIDE)
  
  X_test1 <- X_test
  X_test1$Age <- rep(mean(X_test$Age),nrow(X_test))
  X_test1$SexMale <- rep(mean(X_test$SexMale),nrow(X_test))
  pred <- predict(fit, X_test1)
  perf_no_agesex[i + 2] <- R2(pred = pred, obs = Y_test$CAIDE)
  
}



plotDF <- data.frame(R2 = c(perf_all, perf_agesex, perf_no_agesex),
                     Features = c(rep("All",5),rep("Sex and Age Only",5),rep("No Sex and Age",5)),
                     Method = rep(c("CAIDE1 weights\n(Continuous)", "CAIDE1 weights\n(Discrete)", 
                                    "ElasticNet",
                                    "sPLS", 
                                    "Random Forest"),3)
                     )
plotDF$Method <- factor(plotDF$Method, 
                        levels = rev(c("CAIDE1 weights\n(Discrete)", "CAIDE1 weights\n(Continuous)",
                                       "ElasticNet","sPLS", "Random Forest")))

plotDF$Features <- factor(plotDF$Features,
                          levels = rev(c("All", "Sex and Age Only", "No Sex and Age")))



colors <- c(RColorBrewer::brewer.pal(n = 8, name = "Reds")[4:8],
            RColorBrewer::brewer.pal(n = 8, name = "Oranges")[4:8],
            RColorBrewer::brewer.pal(n = 8, name = "PuRd")[4:8])

# Age
Age <- rep(3, nrow(Y_test))
Age[Y_test$Age < 47] <- 0
Age[Y_test$Age > 53] <- 4

# Sex
Sex <- ifelse(bestModel_test$SexMale$ObsPred_test$predClass == 2,0,1)


maximum <- R2(pred = Sex + Age, obs = Y_test$CAIDE)



p <- ggplot(plotDF) +
  geom_hline(aes(yintercept = maximum), linetype = "dashed", linewidth = 1.5) +
  geom_bar(aes(x = Method, y = R2, fill = Features), 
           stat = "identity", position = position_dodge(), color = "black") +
  coord_flip() +
  theme_bw() +
  xlab(NULL) +
  ylab(expression(R^2)) +
  scale_fill_manual(values = colors) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")


ggsave(p, file = "InfluenceOfSexAge.png", width = 8, height = 5)