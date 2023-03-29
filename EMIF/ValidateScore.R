library(tidyverse)
library(caret)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/metaData_EMIF.RData")
load("EMIF/predictedScore_factors_EMIF.RData")

# Filter data
metaData_EMIF <- as.data.frame(metaData_EMIF)
rownames(metaData_EMIF) <- metaData_EMIF$X
metaData_EMIF <- metaData_EMIF[metaData_EMIF$sex.match == 1,] # Remove unmatching sex
metaData_EMIF <- metaData_EMIF[metaData_EMIF$Age <= 75,] # Midlife

convert <- unique(metaData_EMIF$X[(metaData_EMIF$LastFU_Diagnosis == "MCI") | (metaData_EMIF$LastFU_Diagnosis == "AD")])[-1]
convert1 <- intersect(convert, metaData_EMIF$X[metaData_EMIF$Diagnosis == "NL"])
sci <- metaData_EMIF$X[metaData_EMIF$Diagnosis == "SCI"]
convert1 <- unique(sci, convert1)
metaData_EMIF <- metaData_EMIF[setdiff(rownames(metaData_EMIF),convert1),] # Remove converters and SCI


ABcontrol <- unique(intersect(metaData_EMIF$X[(metaData_EMIF$AMYLOIDstatus == 1)],
                       metaData_EMIF$X[(metaData_EMIF$Diagnosis == "NL")]))
metaData_EMIF <- metaData_EMIF[setdiff(rownames(metaData_EMIF),ABcontrol),] # Remove converters


samples <- intersect(metaData_EMIF$X, rownames(predictedScore_factors))
metaData_fil <- metaData_EMIF[samples,]
save(metaData_fil, file = "EMIF/metaData_fil.RData")
predictedScore_factors_fil <- predictedScore_factors[samples,]
save(predictedScore_factors_fil, file = "EMIF/predictedScore_factors_fil.RData")

load("EMIF/Fit_CombineFactors_CAIDE1_RF.RData")
pred_CAIDE1 <- predict(fit, predictedScore_factors_fil)
load("EMIF/Fit_CombineFactors_LIBRA_RF.RData")
pred_LIBRA <- predict(fit, predictedScore_factors_fil)


testDF <- data.frame(Diagnosis = metaData_fil$Diagnosis,
                     EpiCAIDE1 = pred_CAIDE1,
                     EpiLIBRA = pred_LIBRA,
                     EpiAge = predictedScore_factors_fil$Age,
                     Age = metaData_fil$Age,
                     Sex = metaData_fil$Gender)

testDF$Diagnosis <- ifelse(testDF$Diagnosis == "NL","Control","Cognitive Impaired")
table(testDF$Diagnosis)


score <- c("EpiCAIDE1", "EpiLIBRA", "EpiAge", "Age")
scoreName <- c("Epi-CAIDE1", "Epi-LIBRA", "Epi-Age", "Age")
testDF$Class1 <- factor(testDF$Diagnosis,
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
                      Y = rev(seq(0.1,0.3,length.out = length(aucValue))))

ROCplot$Class <- factor(ROCplot$Class, levels = scoreName)

colors <- rev(c("#E6AB02","#6BAED6","#2171B5","#084594"))
p <- ggplot(ROCplot) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Class), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Score),
            fontface = "bold") +
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

ggsave(p, file = "EMIF/ROC_MCI_EMIF.png", width = 8, height = 5)









model <- lm(Diagnosis ~ Age + EpiAge, data = testDF)
summary(model)
boxplot(pred_CAIDE1[metaData_fil$Diagnosis != "AD"] ~ metaData_fil$MCI[metaData_fil$Diagnosis != "AD"])
boxplot(metaData_fil$Age[metaData_fil$Diagnosis != "AD"] ~ metaData_fil$MCI[metaData_fil$Diagnosis != "AD"])
boxplot(predictedScore_factors_fil$Age[metaData_fil$Diagnosis != "AD"] ~ metaData_fil$MCI[metaData_fil$Diagnosis != "AD"])

t.test(pred_CAIDE1[metaData_fil$Diagnosis == "MCI"],pred_CAIDE1[metaData_fil$Diagnosis == "NL"])
t.test(metaData_fil$Age[metaData_fil$Diagnosis == "MCI"],metaData_fil$Age[metaData_fil$Diagnosis == "NL"])
t.test(predictedScore_factors_fil$Age[metaData_fil$Diagnosis == "MCI"],predictedScore_factors_fil$Age[metaData_fil$Diagnosis == "NL"])

test <- aov(pred_CAIDE1 ~ metaData_fil$Diagnosis)
posthoc <- TukeyHSD(test)
summary(posthoc)



