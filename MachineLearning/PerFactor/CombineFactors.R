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

################################################################################

# Get best model for each factor

################################################################################

# All factors except for age, smoking, and sex
factors <- c("BMI", "Diabetes", "Alcohol", "HDL", "TotalChol", "Physical", "HeartDisease",
             "Education", "Depression", "SysBP", "Diet")
factorNames <- c("BMI", "Type II Diabetes", "L-M Alchol Intake", "HDL Chol.",
                 "Total Chol.", "Physical Act.", "Heart Disease", "Education",
                 "Depression", "Systolic BP", "Dietary Intake")

bestModel_test <- list()
bestModel_CV <- list()
for (i in 1:length(factors)){
  # Literature-based feature selection + EN
  load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
  auc_CV_lit <- finalOutput[[factors[i]]]$AUC
  
  # Correlation-based feature selection + EN
  load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
  auc_CV_EN <- finalOutput[[factors[i]]]$AUC
  
  # Correlation-based feature selection + RF
  load("~/PerFactor/Random Forest/finalOutput_CV.RData")
  auc_CV_RF <- finalOutput[[factors[i]]]$AUC 
  
  best <- which.max(c(auc_CV_lit,auc_CV_EN, auc_CV_RF))
  if (best == 1){
    load("~/PerFactor/Literature/finalOutput_lit_test.RData")
    load("~/PerFactor/Literature/finalOutput_lit_CV.RData")
  }
  if (best == 2){
    load("~/PerFactor/ElasticNet/finalOutput_test.RData")
    load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
  }
  if (best == 3){
    load("~/PerFactor/Random Forest/finalOutput_test.RData")
    load("~/PerFactor/Random Forest/finalOutput_CV.RData")
  }
  bestModel_test[[i]] <- finalOutput_test[[factors[i]]]
  bestModel_CV[[i]] <- finalOutput[[factors[i]]]
}

# Sex
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
bestModel_test[[12]] <- finalOutput_test[["SexMale"]]
load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
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
load("~/PerFactor/ElasticNet/finalOutput_test.RData")
bestModel_test[[15]] <- finalOutput_test$Smoking

load("~/PerFactor/ElasticNet/finalOutput_CV.RData")
bestModel_CV[[15]] <- finalOutput$Smoking


names(bestModel_test) <- c(factors, "SexMale","Age47", "Age53", "Smoking")
names(bestModel_CV) <- c(factors, "SexMale","Age47", "Age53", "Smoking")

################################################################################

# CAIDE1

################################################################################

#******************************************************************************#
# Continuous
#******************************************************************************#
load("~/Data/Y_test.RData")

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

# Predicted score
predCAIDE1 <- Age + BMI + Education + Physical + Sex + SysBP + TotalChol


plotDF <- data.frame(Predicted = predCAIDE1,
                     Observed = Y_test$CAIDE)
p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(y = Observed, x = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("CAIDE1") +
  ylab("Observed Score") +
  xlab("Predicted Score") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(p, file = "ObsVsPred_test_FactorCon_CAIDE1.png", width = 8, height = 4)

RMSE(pred = predCAIDE1, obs = Y_test$CAIDE)
R2(pred = predCAIDE1, obs = Y_test$CAIDE)

#******************************************************************************#
# Categorical
#******************************************************************************#
load("~/Data/Y_test.RData")

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

predCAIDE1 <- Age + BMI + Education + Physical + Sex + SysBP + TotalChol

plotDF <- data.frame(Predicted = predCAIDE1,
                     Observed = Y_test$CAIDE)
p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_jitter(aes(y = Observed, x = Predicted, color = Observed-Predicted), 
              width = 0.05, height = 0.1, size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("CAIDE1") +
  ylab("Observed Score") +
  xlab("Predicted Score") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(p, file = "ObsVsPred_test_FactorCat_CAIDE1.png", width = 8, height = 4)


RMSE(pred = predCAIDE1, obs = Y_test$CAIDE)
R2(pred = predCAIDE1, obs = Y_test$CAIDE)

#******************************************************************************#
# linear model
#******************************************************************************#
load("~/Data/Y_test.RData")

# Age
Age <- bestModel_test$Age47$ObsPred_test$predictedAge

# Sex
Sex <- ifelse(bestModel_test$SexMale$ObsPred_test$predClass == "Yes",1,0)

# BMI
BMI <- bestModel_test$BMI$ObsPred_test$pred

# Education
Education <- bestModel_test$Education$ObsPred_test$pred

# Physical Activity
Physical <- bestModel_test$Physical$ObsPred_test$pred

# Systolic blood pressure
SysBP <- bestModel_test$SysBP$ObsPred_test$pred

# Total Cholesterol
TotalChol <- bestModel_test$TotalChol$ObsPred_test$pred

# Observed
Observed <- Y_test$CAIDE


X <- data.frame(Observed, Age, Sex, BMI, Education, Physical, SysBP, TotalChol)

model <- lm(Observed ~ Age + BMI,data = X)
summary(model)

predCAIDE1 <-predict(model)

plotDF <- data.frame(Predicted = predCAIDE1,
                     Observed = Y_test$CAIDE)
p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(y = Observed, x = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("CAIDE1") +
  ylab("Observed Score") +
  xlab("Predicted Score") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(p, file = "ObsVsPred_test_FactorFit_CAIDE1.png", width = 8, height = 4)


RMSE(pred = predCAIDE1, obs = Y_test$CAIDE)
R2(pred = predCAIDE1, obs = Y_test$CAIDE)


################################################################################

# LIBRA

################################################################################

#******************************************************************************#
# Continuous
#******************************************************************************#
load("~/Data/Y_test.RData")

# Diet
threshold <- bestModel_test$Diet$threshold
power <- log(0.5)/log(threshold)
Diet <- (1-bestModel_test$Diet$ObsPred_test$pred^power)*-1.7

# Physical Activity
threshold <- bestModel_test$Physical$threshold
power <- log(0.5)/log(threshold)
Physical <- (1-bestModel_test$Physical$ObsPred_test$pred^power)*1.1

# Smoking
threshold <- bestModel_test$Smoking$threshold
power <- log(0.5)/log(threshold)
Smoking <- (1-bestModel_test$Smoking$ObsPred_test$pred^power)*1.5

# Alcohol
threshold <- bestModel_test$Alcohol$threshold
power <- log(0.5)/log(threshold)
Alcohol <- (1-bestModel_test$Alcohol$ObsPred_test$pred^power)*-1

# BMI
threshold <- bestModel_test$BMI$threshold
power <- log(0.5)/log(threshold)
BMI <- (1-bestModel_test$BMI$ObsPred_test$pred^power)*1.6

# Depression
threshold <- bestModel_test$Depression$threshold
power <- log(0.5)/log(threshold)
Depression <- (1-bestModel_test$Depression$ObsPred_test$pred^power)*2.1

# Diabetes
threshold <- bestModel_test$Diabetes$threshold
power <- log(0.5)/log(threshold)
Diabetes <- (1-bestModel_test$Diabetes$ObsPred_test$pred^power)*1.3

# Systolic blood pressure
threshold <- bestModel_test$SysBP$threshold
power <- log(0.5)/log(threshold)
SysBP <- (1-bestModel_test$SysBP$ObsPred_test$pred^power)*1.6

# HDL
threshold <- bestModel_test$HDL$threshold
power <- log(0.5)/log(threshold)
HDL <- (1-bestModel_test$HDL$ObsPred_test$pred^power)*1.4

# Heart disease
threshold <- bestModel_test$HeartDisease$threshold
power <- log(0.5)/log(threshold)
HeartDisease <- (1-bestModel_test$HeartDisease$ObsPred_test$pred^power)*1

# Predicted score
predLIBRA <- Diet + Physical + Smoking + Alcohol + BMI + Depression + Diabetes +
  SysBP + HDL + HeartDisease


plotDF <- data.frame(Predicted = predLIBRA,
                     Observed = Y_test$LIBRA)
p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(y = Observed, x = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("LIBRA") +
  ylab("Observed Score") +
  xlab("Predicted Score") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(p, file = "ObsVsPred_test_FactorCon_LIBRA.png", width = 8, height = 4)

RMSE(pred = predLIBRA, obs = Y_test$LIBRA)
R2(pred = predLIBRA, obs = Y_test$LIBRA)

#******************************************************************************#
# Categorical
#******************************************************************************#
load("~/Data/Y_test.RData")

# Diet
Diet <- ifelse(bestModel_test$Diet$ObsPred_test$predClass == "Yes",-1.7,0)

# Physical Activity
Physical <- ifelse(bestModel_test$Physical$ObsPred_test$predClass == "Yes",1.1,0)

# Smoking
Smoking <- ifelse(bestModel_test$Smoking$ObsPred_test$predClass == "Yes",1.5,0)

# Alcohol
Alcohol <- ifelse(bestModel_test$Alcohol$ObsPred_test$predClass == "Yes",-1,0)

# BMI
BMI <- ifelse(bestModel_test$BMI$ObsPred_test$predClass == "Yes",1.6,0)

# Depression
Diet <- ifelse(bestModel_test$Depression$ObsPred_test$predClass == "Yes",2.1,0)

# Diabetes
Diabetes <- ifelse(bestModel_test$Diabetes$ObsPred_test$predClass == "Yes",1.3,0)

# Systolic blood pressure
SysBP <- ifelse(bestModel_test$SysBP$ObsPred_test$predClass == "Yes",1.6,0)

# HDL
HDL <- ifelse(bestModel_test$HDL$ObsPred_test$predClass == "Yes",1.4,0)

# Heart disease
HeartDisease <- ifelse(bestModel_test$HeartDisease$ObsPred_test$predClass == "Yes",1,0)

threshold <- bestModel_test$HeartDisease$threshold
power <- log(0.5)/log(threshold)
HeartDisease <- (1-bestModel_test$HeartDisease$ObsPred_test$pred^power)*1

# Predicted score
predLIBRA <- Diet + Physical + Smoking + Alcohol + BMI + Depression + Diabetes +
  SysBP + HDL + HeartDisease


plotDF <- data.frame(Predicted = predLIBRA,
                     Observed = Y_test$LIBRA)
p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_jitter(aes(y = Observed, x = Predicted, color = Observed-Predicted), 
              width = 0.05, height = 0.1, size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("LIBRA") +
  ylab("Observed Score") +
  xlab("Predicted Score") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(p, file = "ObsVsPred_test_FactorCat_LIBRA.png", width = 8, height = 4)


RMSE(pred = predLIBRA, obs = Y_test$LIBRA)
R2(pred = predLIBRA, obs = Y_test$LIBRA)

#******************************************************************************#
# linear model
#******************************************************************************#
load("~/Data/Y_test.RData")

# Diet
Diet <- bestModel_test$Diet$ObsPred_test$pred

# Physical Activity
Physical <-bestModel_test$Physical$ObsPred_test$pred

# Smoking
Smoking <- bestModel_test$Smoking$ObsPred_test$pred

# Alcohol
Alcohol <- bestModel_test$Alcohol$ObsPred_test$pred

# BMI
BMI <- bestModel_test$BMI$ObsPred_test$pred

# Depression
Depression <- bestModel_test$Depression$ObsPred_test$pred

# Diabetes
Diabetes <- bestModel_test$Diabetes$ObsPred_test$pred

# Systolic blood pressure
SysBP <- bestModel_test$SysBP$ObsPred_test$pred

# HDL
HDL <- bestModel_test$HDL$ObsPred_test$pred

# Heart disease
HeartDisease <-bestModel_test$HeartDisease$ObsPred_test$pred


# Observed
Observed <- Y_test$LIBRA


X <- data.frame(Observed, Diet, Physical, Smoking, Alcohol, BMI, Depression, Diabetes,
                  SysBP, HDL, HeartDisease)

X <- data.frame(Observed, Smoking, Depression)

model <- lm(Observed ~ .,data = X)
summary(model)

predLIBRA <-predict(model)

plotDF <- data.frame(Predicted = predLIBRA,
                     Observed = Y_test$LIBRA)
p <- ggplot(plotDF) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(y = Observed, x = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ggtitle("LIBRA") +
  ylab("Observed Score") +
  xlab("Predicted Score") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

ggsave(p, file = "ObsVsPred_test_FactorFit_LIBRA.png", width = 8, height = 4)


RMSE(pred = predLIBRA, obs = Y_test$LIBRA)
R2(pred = predLIBRA, obs = Y_test$LIBRA)
