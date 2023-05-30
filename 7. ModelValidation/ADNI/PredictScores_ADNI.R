library(tidyverse)
library(caret)
library(ggpubr)

# Clear workspace and console
rm(list = ls())
cat("\014") 

################################################################################

# Predict risk scores in ADNI

################################################################################

# Load data
load("ADNI/X_ADNI_imp.RData")
load("ADNI/MetaData_ADNI.RData")

# Keep midlife samples only:
X_ADNI_fil <- t(X_ADNI_imp)
#keepsamples <- MetaData_baseline$Basename[MetaData_baseline$Age <= 75]
#X_ADNI_fil <- t(X_ADNI_imp)[,keepsamples]

# Calculate scores
Scores <- c("CAIDE1", "CAIDE2", "LIBRA")
Models <- c("EN", "sPLS","RF")
pred_list <- list()


for (i in 1:length(Scores)){
  predictedScore <- matrix(NA, nrow = ncol(X_ADNI_fil), ncol = length(Models))
  colnames(predictedScore) <- Models
  rownames(predictedScore) <- colnames(X_ADNI_fil)
  for (m in 1:length(Models)){
    # Load model
    load(paste0(Scores[i],"_Cor/CV_",Scores[i],"_Cor_", Models[m],".RData"))
    
    # Get features needed for model fitting
    features <- colnames(finalModel$trainingData)[-10001]
    
    # Format ADNI data for model testing
    X_ADNI_test <- X_ADNI_fil[features,]
    sum(is.na(X_ADNI_test))
    
    # Make prediction
    predictedScore[,m] <- predict(finalModel, t(X_ADNI_test))
    
  }
  predictedScore <- as.data.frame(predictedScore)
  pred_list[[i]] <- predictedScore
}

names(pred_list) <- Scores
save(pred_list, file = "ADNI/predictedScores_ADNI.RData")


# Format meta data
MetaData_test <- as.data.frame(MetaData_baseline)
rownames(MetaData_test) <- MetaData_baseline$Basename
MetaData_test <- MetaData_test[colnames(X_ADNI_fil),]



#==============================================================================#
# Correlation with other factors at baseline
#==============================================================================#

# Load predicted scores
load("ADNI/predictedScores_ADNI.RData")

# Format meta data
MetaData_test <- as.data.frame(MetaData_baseline)
rownames(MetaData_test) <- MetaData_baseline$Basename
MetaData_test <- MetaData_test[rownames(pred_list$CAIDE1),]

midlife_samples <- MetaData_test$Basename[MetaData_test$Age <= 75]
MetaData_test <- MetaData_test[midlife_samples,]

pred <- cbind.data.frame(pred_list$CAIDE1[midlife_samples,],
                         pred_list$CAIDE2[midlife_samples,],
                         pred_list$LIBRA[midlife_samples,],
                         MetaData_test$Age,
                         ifelse(MetaData_test$Sex == "M",1,0))

colnames(pred) <- c("CAIDE1_EN","CAIDE1_sPLS","CAIDE1_RF",
                    "CAIDE2_EN","CAIDE2_sPLS","CAIDE2_RF",
                    "LIBRA_EN","LIBRA_sPLS","LIBRA_RF",
                    "Age", "Sex")

pred <- pred[midlife_samples,]


factorNames <- c("APOE4", "FDG", "AV45", "TAU",
             "PTAU", "CDRSB", "ADAS11", "ADAS13", "ADASQ4", "MMSE", 
             "RAVLT.immediate","RAVLT.learning","RAVLT.forgetting","RAVLT.perc.forgetting",
             "LDELTOTAL", "TRABSCOR", "FAQ", "MOCA", "EcogPtMem", "EcogPtLang",
             "EcogPtVisspat", "EcogPtPlan", "EcogPtOrgan", "EcogPtDivatt", "EcogPtTotal",
             "EcogSPMem", "EcogSPLang", "EcogSPVisspat", "EcogSPPlan", "EcogSPOrgan", 
             "EcogSPDivatt", "EcogSPTotal", "Ventricles", "Hippocampus", "WholeBrain",
             "Entorhinal", "Fusiform", "MidTemp", "ICV", "mPACCdigit", "mPACCtrailsB")

factorNames <- c("APOE4", "FDG", "AV45", "TAU",
                 "PTAU", "CDRSB", "ADAS11", "ADAS13", "ADASQ4", "MMSE", 
                 "RAVLT.immediate","RAVLT.learning","RAVLT.forgetting","RAVLT.perc.forgetting",
                 "LDELTOTAL", "TRABSCOR", "FAQ", "MOCA", "EcogPtMem", "EcogPtLang",
                 "EcogPtVisspat", "EcogPtPlan", "EcogPtOrgan", "EcogPtDivatt", "EcogPtTotal",
                 "Ventricles", "Hippocampus", "WholeBrain",
                 "Entorhinal", "Fusiform", "MidTemp", "ICV", "mPACCdigit", "mPACCtrailsB")

factorsDF <- as.data.frame(MetaData_test[,factorNames])

colnames(factorsDF) <- c("APOE4", "FDG", "AV45", "TAU",
                         "PTAU", "CDRSB", "ADAS11", "ADAS13", "ADASQ4", "MMSE", 
                         "RAVLT (immediate)","RAVLT (learning)","RAVLT (forgetting)","RAVLT (% forgetting)",
                         "LDELTOTAL", "TRABSCOR", "FAQ", "MOCA", "EcogPtMem", "EcogPtLang",
                         "EcogPtVisspat", "EcogPtPlan", "EcogPtOrgan", "EcogPtDivatt", "EcogPtTotal",
                         "Ventricles", "Hippocampus", "WholeBrain",
                         "Entorhinal", "Fusiform", "MidTemp", "ICV", "mPACCdigit", "mPACCtrailsB")


correlation <- matrix(NA, nrow = ncol(pred), ncol = ncol(factorsDF))
significance <- matrix(NA, nrow = ncol(pred), ncol = ncol(factorsDF))
for (i in 1:length(factorNames)){
  correlation[,i] <- apply(pred,2, function(x) {cor(x, factorsDF[,i], method = "spearman", 
                        use = "pairwise.complete.obs")})
  significance[,i] <- apply(pred,2, function(x) {cor.test(x, factorsDF[,i], method = "spearman", 
                        use = "pairwise.complete.obs")$p.value})
  
}

colnames(correlation) <- colnames(factorsDF)
rownames(correlation) <- colnames(pred)
colnames(significance) <- colnames(factorsDF)
rownames(significance) <- colnames(pred)

plotDF_cor <- gather(as.data.frame(correlation))
plotDF_cor$key2 <- rep(rownames(correlation),ncol(correlation))
plotDF_sig <- gather(as.data.frame(significance))
plotDF_all <- cbind.data.frame(plotDF_cor, plotDF_sig$value)
plotDF_all$Sig <- ifelse(p.adjust(plotDF_all$`plotDF_sig$value`, "fdr") < 0.05, "Yes", "No")


plotDF_all$Group <- rep("Age/Sex", nrow(plotDF_all))
plotDF_all$Group[str_detect(plotDF_all$key2, "CAIDE1")] <- "Epi-CAIDE1"
plotDF_all$Group[str_detect(plotDF_all$key2, "CAIDE2")] <- "Epi-CAIDE2"
plotDF_all$Group[str_detect(plotDF_all$key2, "LIBRA")] <- "Epi-LIBRA"
plotDF_all$Group <- factor(plotDF_all$Group, levels = c("Epi-CAIDE1", "Epi-CAIDE2", "Epi-LIBRA",
                                                        "Age/Sex"))

plotDF_all$Model <- rep(NA, nrow(plotDF_all))
plotDF_all$Model[str_detect(plotDF_all$key2, "Age")] <- "Age"
plotDF_all$Model[str_detect(plotDF_all$key2, "Sex")] <- "Sex"
plotDF_all$Model[str_detect(plotDF_all$key2, "RF")] <- "RF"
plotDF_all$Model[str_detect(plotDF_all$key2, "sPLS")] <- "sPLS"
plotDF_all$Model[str_detect(plotDF_all$key2, "EN")] <- "EN"
plotDF_all$Model <- factor(plotDF_all$Model, levels = c("EN", "sPLS", "RF", "Age", "Sex"))

plotDF_all$key <- factor(plotDF_all$key, levels = c("APOE4", "FDG", "AV45", "TAU",
                                                    "PTAU", "CDRSB", "ADAS11", "ADAS13", "ADASQ4", "MMSE", 
                                                    "RAVLT (immediate)","RAVLT (learning)","RAVLT (forgetting)","RAVLT (% forgetting)",
                                                    "LDELTOTAL", "TRABSCOR", "FAQ", "MOCA", "EcogPtMem", "EcogPtLang",
                                                    "EcogPtVisspat", "EcogPtPlan", "EcogPtOrgan", "EcogPtDivatt", "EcogPtTotal",
                                                    "Ventricles", "Hippocampus", "WholeBrain",
                                                    "Entorhinal", "Fusiform", "MidTemp", "ICV", "mPACCdigit", "mPACCtrailsB"))

p <- ggplot(plotDF_all) +
  geom_tile(aes(x = key, y = Model, fill = value, color = Sig),
           stat = "identity", width = 0.9, height = 0.9, size = 0.5) +
  facet_grid(cols = vars(Group), scale = "free", space = "free") +
  xlab("") +
  ylab("Machine Learning Model") +
  labs(fill  = "Spearman\nCorrelation")+
  ggtitle("ADNI", subtitle = "Correlations with cognitive biomarkers") +
  coord_flip() +
  theme_bw()+
  scale_color_manual(values = c("grey","black")) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) +
  guides(color = "none")


ggsave(p, file = "ADNI_Correlations_baseline.png", width = 8, height = 8)


#==============================================================================#
# Correlation with other factors at baseline
#==============================================================================#

load("~/ADNI/ADNI_metadata_classes_640.Rdata")
load("ADNI/predictedScores_ADNI.RData")
load("ADNI/MetaData_ADNI.RData")

# Format meta data
MetaData_test <- as.data.frame(MetaData_baseline)
rownames(MetaData_test) <- MetaData_baseline$Basename
MetaData_test <- MetaData_test[rownames(pred_list$CAIDE1),]

midlife_samples <- MetaData_test$Basename[MetaData_test$Age <= 75]
MetaData_test <- MetaData_test[midlife_samples,]

# Get predictions
Predictions <- pred_list$CAIDE1[midlife_samples,]

all(MetaData_test$Basename == rownames(Predictions))
testDF <- data.frame(RID = MetaData_test$RID,
                     Age = MetaData_test$Age,
                     DX = MetaData_test$DX.bl,
                     Predictions)

testDF <- testDF %>%
  group_by(RID) %>%
  reframe(EN = mean(EN),
          sPLS = mean(sPLS),
          RF = mean(RF),
          Age = Age,
          DX = DX,
          RID = RID)

testDF <- unique(testDF)

boxplot(testDF$EN ~ testDF$DX)


plotDF <- inner_join(testDF, ADNI_model_res, by = c("RID" = "RID"))

boxplot(plotDF$RF ~ plotDF$ThreeClass)
boxplot(plotDF$RF ~ ifelse(plotDF$diagbl == "CN", "Control", "Case"))

plot(plotDF$EN, plotDF$Slope)


test <- aov(RF ~ThreeClass, plotDF)
t.test(plotDF$RF ~ ifelse(plotDF$diagbl == "CN", "Control", "Case"), var = TRUE)
var.test(plotDF$RF ~ ifelse(plotDF$diagbl == "CN", "Control", "Case"))
ks.test(plotDF$RF, "pnorm", mean=mean(plotDF$RF), sd=sd(plotDF$RF))

library(ggpubr)
ggqqplot(plotDF$RF)

t.test(plotDF$Age ~ ifelse(plotDF$diagbl == "CN", "Control", "Case"), var = TRUE)


plotDF <- plotDF[plotDF$diagbl != "SMC",]
plotDF$CognitiveImpairment <- ifelse(plotDF$diagbl == "CN", "Control", "Mild Cognitive Impairment")


p <- ggplot(plotDF) +
  geom_boxplot(aes(x = CognitiveImpairment, y = RF, fill = CognitiveImpairment), 
               outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(x = CognitiveImpairment, y = RF, color = CognitiveImpairment), width = 0.1) +
  xlab(NULL) +
  ylab("Epi-CAIDE1") +
  scale_color_manual(values = rev(c("#E41A1C","#377EB8"))) +
  scale_fill_manual(values = rev(c("#E41A1C","#377EB8"))) +
  theme_classic() +
  theme(legend.position = "none")


ggsave(p, file = "MCI_CAIDE1_RF_boxplot.png", width = 6, height = 4)



p <- ggplot(plotDF) +
  geom_boxplot(aes(x = CognitiveImpairment, y = Age, fill = CognitiveImpairment), 
               outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(x = CognitiveImpairment, y = Age, color = CognitiveImpairment), width = 0.1) +
  xlab(NULL) +
  ylab("Age") +
  scale_color_manual(values = rev(c("#E41A1C","#377EB8"))) +
  scale_fill_manual(values = rev(c("#E41A1C","#377EB8"))) +
  theme_classic() +
  theme(legend.position = "none")

ggsave(p, file = "MCI_Age_boxplot.png", width = 6, height = 4)



#==============================================================================#
# Kaplan-Meier (MMSE)
#==============================================================================#

library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
#library(condsurv)
library(tidyverse)

# Load predicted scores
load("ADNI/predictedScores_ADNI.RData")
load("ADNI/MetaData_ADNI.RData")

# Format meta data
MetaData_test <- as.data.frame(MetaData_baseline)
rownames(MetaData_test) <- MetaData_baseline$Basename
MetaData_test <- MetaData_test[rownames(pred_list$CAIDE1),]

midlife_samples <- MetaData_test$Basename[MetaData_test$Age <= 69]
MetaData_test <- MetaData_test[midlife_samples,]

# Format predicted score
dataObj <- data.frame(Basename = rownames(pred_list$CAIDE1),
                      Prediction = pred_list$CAIDE1[,3])

rownames(dataObj) <- dataObj$Basename
dataObj <- dataObj[midlife_samples,]

dataObj <- inner_join(dataObj, MetaData_baseline[, c("Basename", "RID", "Age", "Sex")], by = 
                        c("Basename" = "Basename"))

# Get mean prediction for same sample
dataObj <- dataObj %>%
  group_by(RID) %>%
  reframe(Prediction = mean(Prediction),
          Age = Age,
          Sex = Sex,
          sampleID = RID)

dataObj <- unique(dataObj)


# Add data of other time points
dataObj_all <- inner_join(dataObj, 
                          as.data.frame(MetaData_allTime[, c("RID", "VISCODE","MMSE")]), 
                          by = c("sampleID" = "RID"),
                          multiple = "all")

# Add status: MMSE < 26 -> cognitive decline
dataObj_all$Status <- ifelse(dataObj_all$MMSE < 26, 1, 0)
dataObj_all <- dataObj_all[!is.na(dataObj_all$Status),]

# Format time as continuous variable
dataObj_all$Time <- str_remove(dataObj_all$VISCODE, "m0")
dataObj_all$Time <- str_remove(dataObj_all$Time, "m")
dataObj_all$Time[dataObj_all$Time == "bl"] <- 0
dataObj_all$Time <- as.numeric(dataObj_all$Time)

# Select samples that are measures at time = 36
samples36 <- dataObj_all$RID[dataObj_all$Time ==36]
dataObj_all <- dataObj_all[dataObj_all$RID %in% samples36,]
dataObj_all <- dataObj_all[dataObj_all$Time <= 36,]

# Get cases
dataObj_cases <- dataObj_all[dataObj_all$Status == 1,]
dataObj_cases <- arrange(dataObj_cases, Time)
dataObj_cases <- dataObj_cases[!duplicated(dataObj_cases[, "sampleID"]),]

# Remove cases at D0
dataObj_cases <- dataObj_cases[dataObj_cases$Time != 0,]


# Get controls (non-converters)
dataObj_controls <- dataObj_all[!(dataObj_all$sampleID %in% dataObj_cases$sampleID),]
dataObj_controls <- dataObj_controls[!duplicated(dataObj_controls[, "sampleID"]),]
dataObj_controls$Time <- 36

# Combine cases and controls
dataObj_total <- rbind.data.frame(dataObj_cases, dataObj_controls)

# Split into two categories
t_pred <- quantile(dataObj_total$Prediction,0.5)
dataObj_total$ClassPred <- ifelse(dataObj_total$Prediction >= t_pred, "High Epi-CAIDE1", "Low Epi-CAIDE1")
table(dataObj_total$ClassPred)

t_age <- quantile(dataObj_total$Age,0.5)
dataObj_total$ClassAge <- ifelse(dataObj_total$Age >= t_age, "High Age", "Low Age")
table(dataObj_total$ClassAge)





# Predicted class

dataObj_total$Group <- paste0(dataObj_total$ClassAge, ", ", dataObj_total$ClassPred)
table(dataObj_total$Group)
colors <- c(RColorBrewer::brewer.pal(n = 8, "Reds")[c(5,7)],
            RColorBrewer::brewer.pal(n = 8, "Blues")[c(5,7)])[c(2,1,3,4)]
p <- ggplot(dataObj_total) +
  geom_hline(yintercept = t_pred, linetype = "dashed") +
  geom_vline(xintercept = t_age, linetype = "dashed") +
  geom_point(aes(x = Age, y = Prediction, color = Group)) +
  geom_text(aes(x = 72,y = 4.1, label = paste0("R-Squared = ", 0.43)),
            fontface = 3, size = 3) +
  xlab("Age") +
  ylab("Epi-CAIDE1") +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position =  "right")

ggsave(p, file = "ADNI_CAIDE1vsAge.png", width = 10, height = 6)

caret::R2(dataObj$Age, dataObj$Prediction)


#==============================================================================#
# Predicted score
#==============================================================================#

test <- survfit2(Surv(Time, Status) ~ ClassPred, data = dataObj_total)

plotDF <- data.frame(
  Time = c(0,0,test$time),
  n.risk = c(NA,NA,test$n.risk),
  n.event = c(NA,NA,test$n.event),
  n.censor = c(NA, NA,test$n.censor),
  surv = c(1,1,test$surv),
  Class = factor(c("High Epi-LIBRA", "Low Epi-LIBRA",
                   rep("High Epi-LIBRA",4), rep("Low Epi-LIBRA",4)),
                 levels = c("Low Epi-LIBRA", "High Epi-LIBRA")))


plotDF$Time <- plotDF$Time/12

p <- ggplot(plotDF) +
  geom_step(aes(x = Time, y = surv, group = Class, color = Class),
            linewidth = 1.5, alpha = 0.8) +
  geom_point(aes(x = Time, y = surv, group = Class, color = Class), size = 3,
             shape = 15) +
  geom_text(aes(x = 0.3,y = 0.71, label = paste0("p-value = ", 0.5)),
            fontface = "italic", size = 3) +
  scale_color_manual(values = c("#FB6A4A", "#CB181D")) +
  theme_bw() +
  ylab("Global Cognitive Function\n(MMSE \u2265 26)") +
  xlab("Time (years)") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
  
ggsave(p, file = "KaplanMeier_MMSE_LIBRAScore.png", width = 8, height = 5)
  
# Test for statistical difference
survdiff(Surv(Time, Status) ~ ClassPred, data = dataObj_total)


#==============================================================================#
# Age
#==============================================================================#

test <- survfit2(Surv(Time, Status) ~ ClassAge, data = dataObj_total)

plotDF <- data.frame(
  Time = c(0,0,test$time),
  n.risk = c(NA,NA,test$n.risk),
  n.event = c(NA,NA,test$n.event),
  n.censor = c(NA, NA,test$n.censor),
  surv = c(1,1,test$surv),
  Class = factor(c("High Age", "Low Age",
                   rep("High Age",3), rep("Low Age",4)),
                 levels = c("Low Age", "High Age")))


plotDF$Time <- plotDF$Time/12

p <- ggplot(plotDF) +
  geom_step(aes(x = Time, y = surv, group = Class, color = Class),
            linewidth = 1.5, alpha = 0.8) +
  geom_point(aes(x = Time, y = surv, group = Class, color = Class), size = 3,
             shape = 15) +
  geom_text(aes(x = 0.3,y = 0.71, label = paste0("p-value = ", 0.4)),
            fontface = 3, size = 3) +
  scale_color_manual(values = c("#6BAED6", "#2171B5")) +
  theme_bw() +
  ylab("Global Cognitive Function\n(MMSE \u2265 26)") +
  xlab("Time (years)") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave(p, file = "KaplanMeier_MMSE_Age.png", width = 8, height = 5)

# Test for statistical difference
survdiff(Surv(Time, Status) ~ ClassAge, data = dataObj_total)





p <- survfit2(Surv(Time, Status) ~ ClassPred, data = dataObj_total) %>% 
  ggsurvfit(type = "survival") +
  labs(
    x = "Months",
    y = "Global Cognitive Function"
  ) + 
  add_confidence_interval()




################################################################################

# Adjusted for batch

################################################################################

# Load data
load("ADNI/X_ADNI_adj.RData")
load("ADNI/lumi_dpval_ADNI.RData")
load("ADNI/MetaData_ADNI.RData")

# Keep midlife samples only:
keepsamples <- MetaData_baseline$Basename[MetaData_baseline$Age <= 75]
X_ADNI_fil <- X_ADNI_adj[,keepsamples]

# Load CAIDE1 model (ElasticNet)
load("CAIDE1_Cor/CV_CAIDE1_Cor_EN.RData")

# Get (non-zero) coefficients from model
coefficients <- varImp(finalModel,scale = FALSE)$importance
coefficients <- rownames(coefficients)[coefficients$Overall != 0]

# Keep samples with low detection p-value for these features/probes
test <- lumi_dpval[coefficients,]
keepsamples <- colnames(test)[colSums(test>0.1) == 0]

# Get features needed for model fitting
features <- colnames(finalModel$trainingData)[-10001]

# Format ADNI data for model testing
X_ADNI_test <- X_ADNI_fil[features,colnames(X_ADNI_fil) %in% keepsamples]
sum(is.na(X_ADNI_test))

# Make prediction
pred <- predict(finalModel, t(X_ADNI_test))

# Format meta data
MetaData_test <- as.data.frame(MetaData_baseline)
rownames(MetaData_test) <- MetaData_baseline$Basename
MetaData_test <- MetaData_test[colnames(X_ADNI_test),]



#==============================================================================#
# Correlation with age
#==============================================================================#

load("Data/X_test.RData")
load("Data/Y_test.RData")

# Make prediction
X_test_M <- log2(X_test/(1-X_test))
pred_EXTEND <- predict(finalModel, t(X_test_M[features,]))

plotDF <- data.frame(Pred = c(pred, pred_EXTEND),
                     Age = c(MetaData_test$Age, Y_test$Age),
                     Cohort = c(rep("ADNI", length(pred)), 
                                rep("EXTEND: Test", length(pred_EXTEND))))





# Shading = SE
p <- ggplot(plotDF) +
  geom_point(aes(x = Age, y = Pred), color = "#20262E") +
  geom_smooth(aes(x = Age, y = Pred),
              method='lm', formula= y~x, color = "red", fill = "red") +
  facet_grid(rows = vars(Cohort)) +
  xlab("Age") +
  ylab("Predicted CAIDE1 Score") +
  theme_bw()

ggsave(p, file = "ADNI/AgeVsPred_ADNI_adj.png", width = 8, height = 6)



#==============================================================================#
# Correlation with other factors at baseline
#==============================================================================#


factorNames <- c("Age", "Sex", "APOE4", "FDG", "AV45", "TAU",
                 "PTAU", "CDRSB", "ADAS11", "ADAS13", "ADASQ4", "MMSE", 
                 "RAVLT.immediate","RAVLT.learning","RAVLT.forgetting","RAVLT.perc.forgetting",
                 "LDELTOTAL", "TRABSCOR", "FAQ", "MOCA", "EcogPtMem", "EcogPtLang",
                 "EcogPtVisspat", "EcogPtPlan", "EcogPtOrgan", "EcogPtDivatt", "EcogPtTotal",
                 "EcogSPMem", "EcogSPLang", "EcogSPVisspat", "EcogSPPlan", "EcogSPOrgan", 
                 "EcogSPDivatt", "EcogSPTotal", "Ventricles", "Hippocampus", "WholeBrain",
                 "Entorhinal", "Fusiform", "MidTemp", "ICV", "mPACCdigit", "mPACCtrailsB")

factorsDF <- as.data.frame(MetaData_test[,factorNames])
factorsDF$Sex <- ifelse(factorsDF$Sex == "M", 1, 0)

correlation <- rep(NA, length(factorNames))
significance <- rep(NA, length(factorNames))
for (i in 1:length(factorNames)){
  correlation[i] <- cor(pred, factorsDF[,i], method = "spearman", 
                        use = "pairwise.complete.obs")
  significance[i] <- cor.test(pred, factorsDF[,i], method = "spearman", 
                              use = "pairwise.complete.obs")$p.value
  
}

plotDF_cor <- data.frame(Name = factorNames,
                         Cor = correlation,
                         P = p.adjust(significance, "fdr"),
                         Sig <- ifelse(p.adjust(significance, "fdr") < 0.05, "Yes", "No"))

p <- ggplot(plotDF_cor) +
  geom_bar(aes(x = reorder(Name, Cor), y = Cor, fill = Cor, color = Sig),
           stat = "identity") +
  xlab("") +
  ylab("Spearman Correlation") +
  labs(fill  = "Spearman\nCorrelation")+
  ggtitle("ADNI", subtitle = "Correlations with predicted CAIDE1 score") +
  coord_flip() +
  theme_bw()+
  scale_color_manual(values = c("grey","black")) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1,1))+
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) +
  guides(color = "none")


ggsave(p, file = "ADNI_adj_Correlations_CAIDE1_EN.png", width = 8, height = 8)




#==============================================================================#
# Kaplan-Meier (MMSE)
#==============================================================================#

library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)


# Put Class, diagnosis, time, and ID in a data frame
dataObj <- data.frame(Prediction = pred,
                      sampleID = MetaData_test$RID)

dataObj <- dataObj %>%
  group_by(sampleID) %>%
  reframe(Prediction = mean(Prediction),
          sampleID = sampleID)


# Predicted class
min(dataObj$Prediction)
dataObj$Class <- ifelse(dataObj$Prediction >= 7.5, "High Risk", "Intermediate Risk")
dataObj$Class[dataObj$Prediction < 3.5] <- "Low Risk"

# Get longitudinal data
dataObj <- unique(dataObj)
dataObj_all <- inner_join(dataObj, 
                          as.data.frame(MetaData_allTime[, c("RID", "VISCODE","MMSE")]), 
                          by = c("sampleID" = "RID"),
                          multiple = "all")

dataObj_all$Status <- ifelse(dataObj_all$MMSE < 26, 1, 0)
dataObj_all$Time <- str_remove(dataObj_all$VISCODE, "m0")
dataObj_all$Time <- str_remove(dataObj_all$Time, "m")
dataObj_all$Time[dataObj_all$Time == "bl"] <- 0
dataObj_all$Time <- as.numeric(dataObj_all$Time)

# remove the ones with lower baseline
removeSamples <- dataObj_all$sampleID[(dataObj_all$VISCODE == "bl") &
                                        (dataObj_all$Status != 0)]

dataObj_all <- dataObj_all[!(dataObj_all$sampleID %in% removeSamples),]

dataObj_all <- dataObj_all[dataObj_all$Time <= 72,]

test <- survfit2(Surv(Time, Status) ~ Class, data = dataObj_all)

plotDF <- data.frame(
  Time = test$time,
  n.risk = test$n.risk,
  n.event = test$n.event,
  n.censor = test$n.censor,
  surv = test$surv,
  Class = factor(c(rep("High Risk (8-14)",13), rep("Intermediate Risk (4-7)",13)),
                 levels = c("Intermediate Risk (4-7)", "High Risk (8-14)")))


plotDF$Time <- plotDF$Time/12

p <- ggplot(plotDF) +
  geom_step(aes(x = Time, y = surv, group = Class, color = Class),
            linewidth = 2) +
  scale_color_manual(values = c("#FB6A4A", "#CB181D")) +
  theme_classic() +
  ylab("Global Cognitive Function (MMSE \u2265 26)") +
  xlab("Time (years)") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave(p, file = "KaplanMeier_MMSE_CAIDE1Score.png", width = 8, height = 6)

# Test for statistical difference
survdiff(Surv(Time, Status) ~ Class, data = dataObj_all)






p <- survfit2(Surv(Time, Status) ~ ClassPred, data = dataObj_all) %>% 
  ggsurvfit(type = "survival") +
  labs(
    x = "Months",
    y = "Global Cognitive Function"
  ) + 
  add_confidence_interval()


