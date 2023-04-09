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
                 "Total Chol.", "Physical Inact.", "Heart Disease", "Education",
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

perf <- rep(NA,10)
################################################################################

# CAIDE1

################################################################################

load("~/Data/Y_test.RData")


#******************************************************************************#
# Continuous
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


plotDF_con <- data.frame(Predicted = predCAIDE1,
                     Observed = Y_test$CAIDE)


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

plotDF_cat <- data.frame(Predicted = predCAIDE1,
                     Observed = Y_test$CAIDE)



#******************************************************************************#
# Random forest model
#******************************************************************************#
load("~/Data/Y_test.RData")
load("~/PerFactor/Fit_CombineFactors_CAIDE1_RF.RData")
load("~/PerFactor/predictedScore_factors_EXTEND.RData")

# Observed
Observed <- Y_test$CAIDE

predCAIDE1 <- predict(fit, predictedScore_factors[Y_test$Basename,])

plotDF_rf <- data.frame(Predicted = predCAIDE1,
                     Observed = Y_test$CAIDE)


plotDF_all <- rbind.data.frame(plotDF_cat, plotDF_con, plotDF_rf)
plotDF_all$Group <- c(rep("Discrete",152), rep("Continuous", 152), 
                      rep("Random Forest",152))

p <- ggplot(plotDF_all) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(y = Observed, x = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  facet_grid(rows = vars(Group)) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ylab("Observed Score") +
  xlab("Predicted Score") +
  ggtitle("CAIDE1") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "ObsVsPred_test_factorCom_CAIDE1.png", width = 8, height = 8)

################################################################################

# CAIDE2

################################################################################

load("~/Data/Y_test.RData")



#******************************************************************************#
# Continuous
#******************************************************************************#

# Age
Age47 <- bestModel_test$Age47$ObsPred_test$predictedClass
Age53 <- bestModel_test$Age53$ObsPred_test$predictedClass
Age <- rep(3,length(Age47))
Age[Age47 == "Yes"] <- 0
Age[Age53 == "Yes"] <- 5

# Sex
Sex <- ifelse(bestModel_test$SexMale$ObsPred_test$predClass == "Yes",1,0)

# BMI
threshold <- bestModel_test$BMI$threshold
power <- log(0.5)/log(threshold)
BMI <- (1-bestModel_test$BMI$ObsPred_test$pred^power)*2

# Education
threshold <- bestModel_test$Education$threshold
power <- log(0.5)/log(threshold)
Education <- (1-bestModel_test$Education$ObsPred_test$pred^power)*3

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
TotalChol <- (1-bestModel_test$TotalChol$ObsPred_test$pred^power)*1

# Predicted score
predCAIDE2 <- Age + BMI + Education + Physical + Sex + SysBP + TotalChol


plotDF_con <- data.frame(Predicted = predCAIDE2,
                         Observed = Y_test$CAIDE2)


#******************************************************************************#
# Categorical
#******************************************************************************#

# Age
Age47 <- bestModel_test$Age47$ObsPred_test$predictedClass
Age53 <- bestModel_test$Age53$ObsPred_test$predictedClass
Age <- rep(3,length(Age47))
Age[Age47 == "Yes"] <- 0
Age[Age53 == "Yes"] <- 5

# Sex
Sex <- ifelse(bestModel_test$SexMale$ObsPred_test$predClass == "Yes",1,0)

# BMI
BMI <- ifelse(bestModel_test$BMI$ObsPred_test$predClass == "Yes",2,0)

# Education
Education <- ifelse(bestModel_test$Education$ObsPred_test$predClass == "Yes",3,0)

# Physical activity
Physical <- ifelse(bestModel_test$Physical$ObsPred_test$predClass == "Yes",1,0)

# Systolic Blood pressure
SysBP <- ifelse(bestModel_test$SysBP$ObsPred_test$predClass == "Yes",2,0)

# Total Cholesterol
TotalChol <- ifelse(bestModel_test$TotalChol$ObsPred_test$predClass == "Yes",1,0)

predCAIDE2 <- Age + BMI + Education + Physical + Sex + SysBP + TotalChol


plotDF_cat <- data.frame(Predicted = predCAIDE2,
                         Observed = Y_test$CAIDE2)
#******************************************************************************#
# Random forest model
#******************************************************************************#
load("~/Data/Y_test.RData")
load("~/PerFactor/Fit_CombineFactors_CAIDE2_RF.RData")
load("~/PerFactor/predictedScore_factors_EXTEND.RData")

# Observed
Observed <- Y_test$CAIDE2

predCAIDE2 <- predict(fit, predictedScore_factors[Y_test$Basename,])

plotDF_rf <- data.frame(Predicted = predCAIDE2,
                        Observed = Y_test$CAIDE2)


RMSE(pred = predCAIDE1, obs = Y_test$CAIDE)
R2(pred = predCAIDE1, obs = Y_test$CAIDE)



plotDF_all <- rbind.data.frame(plotDF_cat, plotDF_con, plotDF_rf)
plotDF_all$Group <- c(rep("Discrete",152), rep("Continuous", 152), 
                      rep("Random Forest",152))

p <- ggplot(plotDF_all) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(y = Observed, x = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  facet_grid(rows = vars(Group)) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ylab("Observed Score") +
  xlab("Predicted Score") +
  ggtitle("CAIDE2") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "ObsVsPred_test_factorCom_CAIDE2.png", width = 8, height = 8)
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


plotDF_con <- data.frame(Predicted = predLIBRA,
                     Observed = Y_test$LIBRA)


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


plotDF_cat <- data.frame(Predicted = predLIBRA,
                     Observed = Y_test$LIBRA)

#******************************************************************************#
# Random forest model
#******************************************************************************#
load("~/Data/Y_test.RData")
load("~/PerFactor/Fit_CombineFactors_LIBRA_RF.RData")
load("~/PerFactor/predictedScore_factors_EXTEND.RData")

# Observed
Observed <- Y_test$LIBRA

predLIBRA <- predict(fit, predictedScore_factors[Y_test$Basename,])

plotDF_rf <- data.frame(Predicted = predLIBRA,
                        Observed = Y_test$LIBRA)





plotDF_all <- rbind.data.frame(plotDF_cat, plotDF_con, plotDF_rf)
plotDF_all$Group <- c(rep("Discrete",152), rep("Continuous", 152), 
                      rep("Random Forest",152))

p <- ggplot(plotDF_all) +
  geom_abline(aes(intercept = 0, slope = 1), color = "black", linetype = "dashed", linewidth = 1.5) +
  geom_point(aes(y = Observed, x = Predicted, color = Observed-Predicted), size = 2, alpha = 0.8) +
  facet_grid(rows = vars(Group)) +
  scale_color_gradient2(low = "#000072", mid = "#F49D1A", high = "red", midpoint = 0) +
  ylab("Observed Score") +
  xlab("Predicted Score") +
  ggtitle("LIBRA") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "ObsVsPred_test_factorCom_LIBRA.png", width = 8, height = 8)
