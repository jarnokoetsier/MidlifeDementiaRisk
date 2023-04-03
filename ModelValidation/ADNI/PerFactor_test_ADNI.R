
# Load packages
library(tidyverse)
library(caret)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("ADNI/predictedScore_factors_ADNI.RData")
load("ADNI/MetaData_ADNI.RData")
load("~/ADNI/ADNI_metadata_classes_640.Rdata")


# Select midlife samples
MetaData_test <- as.data.frame(MetaData_baseline)
rownames(MetaData_test) <- MetaData_test$Basename
midlife_samples <- intersect(MetaData_test$Basename[MetaData_test$Age <= 92], 
                             rownames(predictedScore_factors))


# Filter for midlife samples
MetaData_test <- MetaData_test[midlife_samples,]
X_test <-  predictedScore_factors[midlife_samples,]


# Add PGS
PGS <- read.delim("ADNI/PGS_all_ADNI.tsv")
common_samples <- intersect(MetaData_test$PTID, PGS$X)
MetaData_test <- inner_join(MetaData_test, PGS, by = c("PTID" = "X"))
X_test <- X_test[MetaData_test$Basename,]

# Make prediction
load("~/Fit_CombineFactors_CAIDE1_RF.RData")
pred_CAIDE1 <- predict(fit, X_test)
load("~/Fit_CombineFactors_LIBRA_RF.RData")
pred_LIBRA <- predict(fit, X_test)

# Combine into data frame
testDF <- data.frame(RID = MetaData_test$RID,
                     Age = MetaData_test$Age,
                     Sex = MetaData_test$Sex,
                     DX = MetaData_test$DX.bl,
                     EpiCAIDE1 = pred_CAIDE1,
                     EpiLIBRA = pred_LIBRA,
                     AD_noAPOE = MetaData_test$Ad_no_APOE_bayesr,
                     AD1 = MetaData_test$AD1_bayesr,
                     AD2 = MetaData_test$AD2_bayesr)

# Get mean predicted score for the same sample
testDF <- testDF %>%
  group_by(RID) %>%
  reframe(EpiCAIDE1 = mean(EpiCAIDE1),
          EpiLIBRA = mean(EpiLIBRA),
          Age = Age,
          Sex = Sex,
          DX = DX,
          RID = RID,
          AD_noAPOE = AD_noAPOE,
          AD1 = AD1,
          AD2 = AD2)

testDF <- unique(testDF)


# Combine with filtered meta data
plotDF <- inner_join(testDF, ADNI_model_res, by = c("RID" = "RID"))
plotDF <- plotDF[plotDF$diagbl != "SMC", ]
plotDF$MCI <- ifelse(plotDF$diagbl == "CN", "Control", "MCI")
table(plotDF$MCI)



#==============================================================================#
# Correlation at baseline
#==============================================================================#


# Get cognitive biomarkers
factorNames <- c("FDG", "AV45", "TAU",
                 "PTAU", "CDRSB", "ADAS11", "ADAS13", "ADASQ4", "MMSE", 
                 "RAVLT.immediate","RAVLT.learning","RAVLT.forgetting","RAVLT.perc.forgetting",
                 "LDELTOTAL", "TRABSCOR", "FAQ", "MOCA", "EcogPtMem", "EcogPtLang",
                 "EcogPtVisspat", "EcogPtPlan", "EcogPtOrgan", "EcogPtDivatt", "EcogPtTotal",
                 "Ventricles", "Hippocampus", "WholeBrain",
                 "Entorhinal", "Fusiform", "MidTemp", "ICV", "mPACCdigit", "mPACCtrailsB")



# Get samples
factorsDF <- unique(MetaData_test[,c("RID",factorNames)])
plotDF_factors <- inner_join(plotDF[,c("RID")], factorsDF, 
                            by = c("RID" = "RID"))
all(plotDF_factors$RID == plotDF$RID)

colnames(plotDF_factors) <- c("RID","FDG", "AV45", "TAU",
                         "PTAU", "CDRSB", "ADAS11", "ADAS13", "ADASQ4", "MMSE", 
                         "RAVLT (immediate)","RAVLT (learning)","RAVLT (forgetting)","RAVLT (% forgetting)",
                         "LDELTOTAL", "TRABSCOR", "FAQ", "MOCA", "EcogPtMem", "EcogPtLang",
                         "EcogPtVisspat", "EcogPtPlan", "EcogPtOrgan", "EcogPtDivatt", "EcogPtTotal",
                         "Ventricles", "Hippocampus", "WholeBrain",
                         "Entorhinal", "Fusiform", "MidTemp", "ICV", "mPACCdigit", "mPACCtrailsB")

plotDF_factors <- as.data.frame(plotDF_factors)

# Calculate correlations
pred <- as.data.frame(plotDF[, c("EpiCAIDE1", "EpiLIBRA", "AD1","AD_noAPOE","Age", "Sex")])
pred$Sex <- ifelse(pred$Sex == "M",1,0)
colnames(pred) <- c("Epi-CAIDE1", "Epi-LIBRA", "AD",
                    "AD\n(w/o APOE)","Age", "Sex")
  
  
correlation <- matrix(NA, nrow = ncol(pred), ncol = ncol(plotDF_factors)-1)
significance <- matrix(NA, nrow = ncol(pred), ncol = ncol(plotDF_factors)-1)
for (i in 2:ncol(plotDF_factors)){
  correlation[,i-1] <- apply(pred,2, function(x) {cor(x, plotDF_factors[,i], method = "spearman", 
                                                    use = "pairwise.complete.obs")})
  significance[,i-1] <- apply(pred,2, function(x) {cor.test(x, plotDF_factors[,i], method = "spearman", 
                                                          use = "pairwise.complete.obs")$p.value})
  
}

colnames(correlation) <- colnames(plotDF_factors[,-1])
rownames(correlation) <- colnames(pred)
colnames(significance) <- colnames(plotDF_factors[,-1])
rownames(significance) <- colnames(pred)

plotDF_cor <- gather(as.data.frame(correlation))
plotDF_cor$key2 <- rep(rownames(correlation),ncol(correlation))
plotDF_sig <- gather(as.data.frame(significance))
plotDF_all <- cbind.data.frame(plotDF_cor, plotDF_sig$value)
plotDF_all$Sig <- ifelse(p.adjust(plotDF_all$`plotDF_sig$value`, "fdr") < 0.05, "Yes", "No")


plotDF_all$key <- factor(plotDF_all$key, levels = c("FDG", "AV45", "TAU",
                                                    "PTAU", "CDRSB", "ADAS11", "ADAS13", "ADASQ4", "MMSE", 
                                                    "RAVLT (immediate)","RAVLT (learning)","RAVLT (forgetting)","RAVLT (% forgetting)",
                                                    "LDELTOTAL", "TRABSCOR", "FAQ", "MOCA", "EcogPtMem", "EcogPtLang",
                                                    "EcogPtVisspat", "EcogPtPlan", "EcogPtOrgan", "EcogPtDivatt", "EcogPtTotal",
                                                    "Ventricles", "Hippocampus", "WholeBrain",
                                                    "Entorhinal", "Fusiform", "MidTemp", "ICV", "mPACCdigit", "mPACCtrailsB"))

plotDF_all$Group <- rep(c("Epigenetic", "Epigenetic", "Genetic",
                          "Genetic","Demographic", "Demographic"), ncol(correlation))

plotDF_all$Group <- factor(plotDF_all$Group, levels = c("Epigenetic", "Genetic", "Demographic"))

p <- ggplot(plotDF_all) +
  geom_tile(aes(x = key, y = key2, fill = value, color = Sig),
            stat = "identity", width = 0.9, height = 0.9, size = 0.5) +
  facet_grid(cols = vars(Group), scale = "free", space = "free") +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill  = "Spearman\nCorrelation")+
  ggtitle("ADNI", subtitle = "Correlations with cognitive biomarkers") +
  coord_flip() +
  theme_bw()+
  scale_color_manual(values = c("grey","black")) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       oob = scales::squish,
                       limits = c(-0.5,0.5))+
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) +
  guides(color = "none")


ggsave(p, file = "ADNI_Correlations_baseline.png", width = 8, height = 8)




#==============================================================================#
# Prediction of MCI
#==============================================================================#
table(plotDF$MCI)
t.test(plotDF$EpiCAIDE1 ~ plotDF$MCI, var = FALSE)

colors <- RColorBrewer::brewer.pal(n = 8, name = "Set1")[c(2,1)]
plotDF$MCI_name <- ifelse(plotDF$MCI == "MCI", "Mild Cognitive Impairment", "Healthy Control")
p <- ggplot(plotDF) +
  geom_boxplot(aes(x = MCI_name, y = EpiCAIDE1, fill = MCI), 
               outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(x = MCI_name, y = EpiCAIDE1, color = MCI), 
              width = 0.1) +
  xlab(NULL) +
  ylab("Epi-CAIDE1") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(p, file = "Boxplot_MCI_CAIDE1_ADNI.png", width = 8, height = 5)

library(pROC)

score <- c("EpiCAIDE1", "EpiLIBRA", "AD1", "AD_noAPOE", "Age")
scoreName <- c("Epi-CAIDE1", "Epi-LIBRA", "AD", "AD (w/o APOE)", "Age")

plotDF <- as.data.frame(plotDF)
ROCplot <- NULL
aucValue <- rep(NA, length(score))
for (i in 1:length(score)){
  test <- roc(plotDF$MCI, plotDF[,score[i]])
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- round(as.numeric(auc(test)),2)
}

plotAUC <- data.frame(AUC = paste0("AUC = ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.35,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- rev(c("#E6AB02","#FB6A4A","#CB181D","#6BAED6","#2171B5"))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                                      color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  scale_color_manual(values = colors) +
  ggtitle("MCI") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "ROC_MCI_ADNI.png", width = 8, height = 5)
  


#==============================================================================#
# Prediction of Converter
#==============================================================================#

table(plotDF$conversion)
score <- c("EpiCAIDE1", "EpiLIBRA", "AD1", "AD_noAPOE", "Age")
scoreName <- c("Epi-CAIDE1", "Epi-LIBRA", "AD", "AD (w/o APOE)", "Age")

plotDF <- as.data.frame(plotDF)
ROCplot <- NULL
aucValue <- rep(NA, length(score))
for (i in 1:length(score)){
  test <- roc(factor(plotDF$conversion, levels = c("non-converter", "converter")), 
              plotDF[,score[i]])
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- round(as.numeric(auc(test)),2)
}

plotAUC <- data.frame(AUC = paste0("AUC = ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.35,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- rev(c("#E6AB02","#FB6A4A","#CB181D","#6BAED6","#2171B5"))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  scale_color_manual(values = colors) +
  ggtitle("Converter") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "ROC_converter_ADNI.png", width = 8, height = 5)


#==============================================================================#
# Prediction of Decliner
#==============================================================================#

table(plotDF$TwoClass)
score <- c("EpiCAIDE1", "EpiLIBRA", "AD1", "AD_noAPOE", "Age")
scoreName <- c("Epi-CAIDE1", "Epi-LIBRA", "AD", "AD (w/o APOE)", "Age")

plotDF$TwoClass_name <- ifelse(plotDF$TwoClass == 1, "Decliner", "Stable")
plotDF <- as.data.frame(plotDF)
ROCplot <- NULL
aucValue <- rep(NA, length(score))
for (i in 1:length(score)){
  test <- roc(factor(plotDF$TwoClass_name, levels = c("Stable", "Decliner")), 
              plotDF[,score[i]])
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
}

plotAUC <- data.frame(AUC = paste0("AUC = ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.35,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- rev(c("#E6AB02","#FB6A4A","#CB181D","#6BAED6","#2171B5"))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  scale_color_manual(values = colors) +
  ggtitle("Decliner") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "ROC_decliner_ADNI.png", width = 8, height = 5)












plotDF$conversion_num <- ifelse(plotDF$conversion == "converter",1,0)
plotDF$MCI_num <- ifelse(plotDF$Test == "MCI",1,0)
model <- lm(MCI_num ~ Pred + Age, data = plotDF)
summary(model)




table(testDF$DX)
boxplot(testDF$Pred ~ testDF$DX)
test <- aov(Pred ~DX, testDF)
summary(test)

testDF_fil <- testDF[(testDF$DX == "EMCI") | (testDF$DX == "CN") | (testDF$DX == "LMCI"),]
t.test(testDF_fil$Pred ~ ifelse(testDF_fil$DX == "CN", "Control", "Case"), var = TRUE)
boxplot(testDF_fil$Pred ~ testDF_fil$DX)

ADNI_older_fil <- ADNI_older[ADNI_older$VISCODE == "bl",]
testDF$RID <- as.character(testDF$RID)
plotDF <- inner_join(testDF, ADNI_older_fil, by = c("RID" = "RID"))





boxplot(plotDF$Pred ~ plotDF$ThreeClass)
boxplot(plotDF$Pred ~ ifelse(plotDF$DX == "CN", "Control", "Case"))

plot(plotDF$EN, plotDF$Slope)

plotDF <- plotDF[plotDF$diagbl != "SMC",]
test <- aov(Pred ~ThreeClass, plotDF)
t.test(plotDF$Pred ~ ifelse(plotDF$diagbl == "CN", "Control", "Case"), var = FALSE)
var.test(plotDF$Pred ~ ifelse(plotDF$diagbl == "CN", "Control", "Case"))
ks.test(plotDF$Pred, "pnorm", mean=mean(plotDF$Pred), sd=sd(plotDF$Pred))

library(ggpubr)
ggqqplot(plotDF$RF)

t.test(plotDF$Age ~ ifelse(plotDF$diagbl == "CN", "Control", "Case"), var = TRUE)



plotDF$CognitiveImpairment <- ifelse(plotDF$diagbl == "CN", "Control", "Mild Cognitive Impairment")



#*****************************************************************************#
# Baseline predictions
#*****************************************************************************#

# Format meta data
MetaData_test <- as.data.frame(MetaData_baseline)
rownames(MetaData_test) <- MetaData_baseline$Basename
MetaData_test <- MetaData_test[rownames(X_test),]

pred <- as.data.frame(testPred)

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


plotDF_all$key <- factor(plotDF_all$key, levels = c("APOE4", "FDG", "AV45", "TAU",
                                                    "PTAU", "CDRSB", "ADAS11", "ADAS13", "ADASQ4", "MMSE", 
                                                    "RAVLT (immediate)","RAVLT (learning)","RAVLT (forgetting)","RAVLT (% forgetting)",
                                                    "LDELTOTAL", "TRABSCOR", "FAQ", "MOCA", "EcogPtMem", "EcogPtLang",
                                                    "EcogPtVisspat", "EcogPtPlan", "EcogPtOrgan", "EcogPtDivatt", "EcogPtTotal",
                                                    "Ventricles", "Hippocampus", "WholeBrain",
                                                    "Entorhinal", "Fusiform", "MidTemp", "ICV", "mPACCdigit", "mPACCtrailsB"))

p <- ggplot(plotDF_all) +
  geom_tile(aes(x = key, y = key2, fill = value, color = Sig),
            stat = "identity", width = 0.9, height = 0.9, size = 0.5) +
  #facet_grid(cols = vars(Group), scale = "free", space = "free") +
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


ggsave(p, file = "ADNI_Correlations_CAIDE1_EN.png", width = 8, height = 8)



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

# Format meta data
MetaData_test <- as.data.frame(MetaData_baseline)
rownames(MetaData_test) <- MetaData_baseline$Basename
MetaData_test <- MetaData_test[rownames(X_test),]

midlife_samples <- MetaData_test$Basename[MetaData_test$Age < 75]
MetaData_test <- MetaData_test[midlife_samples,]

# Format predicted score
X_test <-  predictedScore_factors[midlife_samples,]
testPred <- predict(fit,X_test)
pred <- as.data.frame(testPred)
dataObj <- data.frame(Basename = rownames(X_test),
                      Prediction = testPred)

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
                   rep("High Age",4), rep("Low Age",4)),
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


