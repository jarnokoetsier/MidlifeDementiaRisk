# Load packages
library(tidyverse)
library(caret)
library(spls)
library(ranger)
library(glmnet)
library(pROC)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("PPMI/predictedScore_factors_PPMI.RData")
load("PPMI/metaData_ppmi.RData")

# Combine scores and meta data
metaData_all <- metaData_all[(metaData_all$age >= 40) & (metaData_all$age <= 75),]
samples <- intersect(metaData_all$Basename, rownames(predictedScore_factors))
rownames(metaData_all) <- metaData_all$Basename
metaData_fil <- metaData_all[samples,]
predictedScore_factors_fil <- predictedScore_factors[samples,]

length(unique(metaData_fil$ID))

# Make prediction
load("~/Fit_CombineFactors_CAIDE1_RF.RData")
pred_CAIDE1 <- predict(fit, predictedScore_factors_fil)
load("~/Fit_CombineFactors_LIBRA_RF.RData")
pred_LIBRA <- predict(fit, predictedScore_factors_fil)
load("~/PPMI/Fit_EMIF_MCI_EN.RData")
pred_EN <- predict(fit, predictedScore_factors_fil, type = "prob")
load("~/PPMI/Fit_EMIF_MCI_sPLS.RData")
pred_sPLS <- predict(fit, predictedScore_factors_fil, type = "prob")
load("~/PPMI/Fit_EMIF_MCI_RF.RData")
pred_RF <- predict(fit, predictedScore_factors_fil, type = "prob")

# Make boxplots
boxplot(pred_RF$MCI ~ as.character(metaData_fil$Class))
boxplot(pred_CAIDE1 ~ as.character(metaData_fil$Class))
boxplot(pred_LIBRA ~ as.character(metaData_fil$Class))
boxplot(predictedScore_factors_fil$Age ~ as.character(metaData_fil$Class))
boxplot(metaData_fil$age ~ as.character(metaData_fil$Class))

# Perform t-tests
t.test(pred_CAIDE1 ~ as.character(metaData_fil$Class))
t.test(pred_LIBRA ~ as.character(metaData_fil$Class))
t.test(predictedScore_factors_fil$Age ~ as.character(metaData_fil$Class))
t.test(metaData_fil$age ~ as.character(metaData_fil$Class))
t.test(pred_RF$MCI ~ as.character(metaData_fil$Class))

# Combine into data frame
testDF <- cbind.data.frame(predictedScore_factors_fil,
                           pred_CAIDE1,
                           pred_LIBRA,
                           pred_EN = pred_EN$MCI,
                           pred_sPLS = pred_sPLS$MCI,
                           pred_RF = pred_RF$MCI,
                           metaData_fil)

# Adjust for age and sex as covariates
model <- lm(Class ~ age + sex  + pred_EN, data = testDF)
summary(model)

#*****************************************************************************#
#   Make ROC curves
#*****************************************************************************#

score <- c("pred_EN", "pred_sPLS","pred_RF","pred_CAIDE1", "pred_LIBRA", "Age", "age")
scoreName <- c("MCI model (EN)","MCI model (sPLS-DA)","MCI model (RF)","Epi-CAIDE1", "Epi-LIBRA", "Epi-Age", "Age")
testDF$Class1 <- factor(ifelse(testDF$Class == 1, "Cognitive Impaired", "Control"),
                        levels = c("Control", "Cognitive Impaired"))
plotDF <- as.data.frame(testDF)
ROCplot <- NULL
aucValue <- rep(NA, length(score))
for (i in 1:length(score)){
  test <- roc(plotDF$Class1, plotDF[,score[i]])
  
  temp <- data.frame(Sensitivity = test$sensitivities,
                     Specificity = test$specificities,
                     Class = rep(scoreName[i],length(test$specificities)))
  
  ROCplot <- rbind.data.frame(ROCplot, temp)
  aucValue[i] <- format(round(as.numeric(auc(test)),2),nsmall = 2)
}

plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Score = scoreName,
                      X = 0.9,
                      Y = rev(seq(0.05,0.4,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

# Set colors
colors <- rev(c("#E6AB02","#6BAED6","#2171B5","#084594","#EF3B2C","#CB181D", "#99000D"))

# Make plot
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
  ggtitle("PPMI: 8-year cognitive outcome") +
  scale_color_manual(values = colors) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "ROC_CI_PPMI_v2.png", width = 8, height = 5)



#*****************************************************************************#
#   Correlations at baseline
#*****************************************************************************#

predDF <- testDF[,c("pred_EN", "pred_sPLS","pred_RF","pred_CAIDE1", "pred_LIBRA", "Age", "age", "sex")]
colnames(predDF) <- c("MCI model (EN)","MCI model (sPLS-DA)","MCI model (RF)",
                      "Epi-CAIDE1", "Epi-LIBRA", "Epi-Age", "Age", "Sex")

obsDF <- testDF[,c("abeta", "tau", "ptau", "tau_ab", "ptau_ab", 
                   "hvlt_immediaterecall", "HVLTRDLY", "hvlt_discrimination",
                   "HVLTREC", "hvlt_retention","bjlot", 
                   "SDMTOTAL", "sft","VLTANIM", "VLTVEG", "VLTFRUIT",
                   "stai", "stai_state", "stai_trait","upsit", "moca")]

colnames(obsDF) <- c("A-Beta", "Tau", "p-Tau", "Tau/A-Beta ratio", "p-Tau/A-Beta ratio",
                     "HVLT (immediate recall)", "HVLT (delayed recall)", "HVLT (discrimination)",
                     "HVLT (recognition)", "HVLT (retention)", "JoLO Test Score",
                     "Symbol Digit Modalities Score", "Semantic Fluency Score (Total)",
                     "Semantic Fluency Score (Animals)","Semantic Fluency Score (Vegatables)",
                     "Semantic Fluency Score (Fruits)", "STAI Total Score", "STAI State Sub-score", 
                     "STAI Trait Sub-score", "UPSIT Score", "MoCA Score")

# Calculate correlations and significance
correlation <- matrix(NA, nrow = ncol(predDF), ncol = ncol(obsDF))
significance <- matrix(NA, nrow = ncol(predDF), ncol = ncol(obsDF))
for (i in 1:ncol(obsDF)){
  correlation[,i] <- apply(predDF,2, function(x) {cor(x, obsDF[,i], method = "spearman", 
                                                      use = "pairwise.complete.obs")})
  significance[,i] <- apply(predDF,2, function(x) {cor.test(x, obsDF[,i], method = "spearman", 
                                                            use = "pairwise.complete.obs")$p.value})
  
}

colnames(correlation) <- colnames(obsDF)
rownames(correlation) <- colnames(predDF)
colnames(significance) <- colnames(obsDF)
rownames(significance) <- colnames(predDF)


# Format data for plotting
plotDF_cor <- gather(as.data.frame(correlation))
plotDF_cor$key2 <- rep(rownames(correlation),ncol(correlation))
plotDF_sig <- gather(as.data.frame(significance))
plotDF_all <- cbind.data.frame(plotDF_cor, plotDF_sig$value)
plotDF_all$Sig <- ifelse(p.adjust(plotDF_all$`plotDF_sig$value`, "fdr") < 0.05, "Yes", "No")


plotDF_all$Group <- rep(c("Epigenetic", "Epigenetic", "Epigenetic", "Epigenetic", "Epigenetic", "Epigenetic",
                          "Demographic", "Demographic"), ncol(correlation))

plotDF_all$Group <- factor(plotDF_all$Group, levels = c("Epigenetic", "Genetic", "Demographic"))

plotDF_all$key <- factor(plotDF_all$key, levels = 
                           c("A-Beta", "Tau", "p-Tau", "Tau/A-Beta ratio", "p-Tau/A-Beta ratio",
                             "HVLT (immediate recall)", "HVLT (delayed recall)", "HVLT (discrimination)",
                             "HVLT (recognition)", "HVLT (retention)", "JoLO Test Score",
                             "Symbol Digit Modalities Score", "Semantic Fluency Score (Total)",
                             "Semantic Fluency Score (Animals)","Semantic Fluency Score (Vegatables)",
                             "Semantic Fluency Score (Fruits)", "STAI Total Score", "STAI State Sub-score", 
                             "STAI Trait Sub-score", "UPSIT Score", "MoCA Score"))

plotDF_all$key2 <- factor(plotDF_all$key2, levels = c("MCI model (EN)","MCI model (sPLS-DA)","MCI model (RF)",
                                                      "Epi-CAIDE1", "Epi-LIBRA", "Epi-Age", "Age", "Sex"))
# Make correlaiton heatmap
p <- ggplot(plotDF_all) +
  geom_tile(aes(x = key, y = key2, fill = value, color = Sig),
            stat = "identity", width = 0.9, height = 0.9, size = 0.5) +
  facet_grid(cols = vars(Group), scale = "free", space = "free") +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill  = "Spearman\nCorrelation")+
  ggtitle("PPMI", subtitle = "Correlations with cognitive biomarkers") +
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
                                     face = "italic"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  guides(color = "none")

# Save plot
ggsave(p, file = "PPMI_Correlations_baseline.png", width = 9, height = 8)

# Other correlations:
cor.test(1-testDF$Depression, testDF$NP1DPRS, method = "spearman")
cor.test(1-testDF$Diet, testDF$quip_eat, method = "spearman")
cor.test(testDF$moca, testDF$Education, method = "spearman")
cor.test(testDF$moca, testDF$EDUCYRS, method = "spearman")
cor.test(testDF$Education, testDF$EDUCYRS, method = "spearman")
cor.test(1-testDF$HeartDisease, testDF$scopa_cv, method = "spearman")
