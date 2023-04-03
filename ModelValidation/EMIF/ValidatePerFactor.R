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
metaData_EMIF <- metaData_EMIF[metaData_EMIF$sex.match == 1,] # Remove unmatching sex
rownames(metaData_EMIF) <- metaData_EMIF$X
samples <- intersect(metaData_EMIF$X, rownames(predictedScore_factors))
metaData_fil <- metaData_EMIF[samples,]
predictedScore_factors_fil <- predictedScore_factors[samples,]


PhenoNames <- c("Hypertension", "Alcohol Consumption", "Cardiovascular Disease",
                "Depression", "Years of Education", "Smoking Status", "Statines")
testDF <- metaData_fil[,c("Hypertension", "Alcohol", "CardiovascularDis", 
                          "Depression", "Eduy", "Smoking",
                          "Statines")]

temp <- predictedScore_factors_fil
temp <- 1-temp
temp$Education <- predictedScore_factors_fil$Education
temp$Age <- predictedScore_factors_fil$Age
predictedScore_factors_fil <- temp
rm(temp)

MRS_test <- setdiff(colnames(predictedScore_factors_fil), c("SexMale", "Age"))
MRS_names <- c("BMI", "Type II Diabetes", "L-M Alcohol Intake", "HDL Cholesterol",
               "Total Cholesterol", "Physical Inactivity", "Heart Disease", "Education",
               "Depression", "Systolic Blood Pressure", "Dietary Intake", "Smoking")
corDF <- NULL
for (i in 1:ncol(testDF)){
  temp_cor <- apply(predictedScore_factors_fil[,MRS_test],2,function(x){cor(x,as.numeric(testDF[,i]),
                                                     method = "spearman", use = "pairwise.complete.obs")})
  temp_p <- apply(predictedScore_factors_fil[,MRS_test],2,function(x){cor.test(x,as.numeric(testDF[,i]),
                                                     method = "spearman", use = "pairwise.complete.obs")$p.value})
  
  temp <- data.frame(MRS = MRS_names,
                     Pheno = rep(PhenoNames[i],length(temp_cor)),
                     Cor = temp_cor,
                     Pvalue = temp_p)
  corDF <- rbind.data.frame(corDF, temp)
}

corDF$FDR <- p.adjust(corDF$Pvalue, method = "fdr")
corDF$Sig <- ifelse(corDF$FDR < 0.05, "Yes", "No")

p <- ggplot(corDF) +
  geom_tile(aes(x = MRS, y = Pheno, fill = Cor, color = Sig),
            stat = "identity", width = 0.9, height = 0.9, size = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill  = "Spearman\nCorrelation")+
  #coord_flip() +
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


ggsave(p, file = "EMIF/EMIF_Correlations_MRS.png", width = 7, height = 6)




library(pROC)

temp <- metaData_fil
temp$Alcohol <- ifelse(temp$Smoking == 0, 0,1)
temp$Smoking <- ifelse(temp$Smoking == 0, 0,1)
temp$Eduy <- ifelse(temp$Eduy <= 9,1,0)
test <- data.frame(Pheno = c("Depression", "Hypertension", "Alcohol", "Smoking", "Eduy"),
                   MRS = c("Depression", "SysBP", "Alcohol", "Smoking", "Education"))

plotDF <- NULL
aucValue <- rep(NA, nrow(test))
name123 <- c("Depression", "Hypertension", "Alcohol Consumption", "Smoking Status", "Education")
for (i in 1:nrow(test)){
  r <- roc(temp[,test$Pheno[i]], predictedScore_factors_fil[, test$MRS[i]])
  temp_roc <- data.frame(Sensitivity = r$sensitivities,
                         Specificity = r$specificities,
                         Name = rep(name123[i],length(r$specificities)))
  aucValue[i] <- format(round(as.numeric(auc(r)),2),nsmall = 2)
  
  plotDF <- rbind.data.frame(plotDF, temp_roc)
}
plotAUC <- data.frame(AUC = paste0("AUC: ",aucValue),
                      Name = name123,
                      X = 0.9,
                      Y = c(0.25,0.2,0.15,0.1,0.05))
plotDF$Name <- factor(plotDF$Name, levels = name123)

p <- ggplot(plotDF) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 2) +
  geom_path(aes(y = Sensitivity, x = 1- Specificity,
                color = Name), 
            linewidth = 1.5, linetype = "solid") +
  geom_text(data = plotAUC, aes(x = X, y = Y, label = AUC, color = Name),
            fontface = "bold") +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

ggsave(p, file = "EMIF/ROC_factors_EMIF.png", width = 8, height = 5)

test <- roc(metaData_fil$Depression, predictedScore_factors_fil$Depression)
