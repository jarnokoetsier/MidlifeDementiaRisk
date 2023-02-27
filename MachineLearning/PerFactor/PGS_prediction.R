
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
load("~/PGS/df_list.RData")
load("~/Data/metaData_ageFil.RData")
load("PerFactor/Y_test_factors.RData")
load("PerFactor/Y_nonTest_factors.RData")

Y_all <- rbind.data.frame(Y_test, Y_nonTest)
Y_all$Basename <- rownames(Y_all)
Y_all <- inner_join(dat[,c("ID", "Basename")], Y_all, by = c("Basename" = "Basename"))
rownames(Y_all) <- Y_all$ID

samples <- intersect(rownames(df_list$lasso), rownames(Y_all))


PGS_factors <- list(
  BMI = c("BMI", "ExtremeBMI", "ObesityClass1", "ObesityClass2", 
          "ObesityClass3", "Overweight"),
  Diabetes = c("T2D"),
  Alcohol = c("DC2", "Alcohol"),
  Depression = c("MDD"),
  Diet = c("DC2"),
  Education = c("EA", "EA22"),
  HDL = c("HDL"),
  HeartDisease = c("CAD"),
  Physical = c("LST", "MVPA"),
  SysBP = c("SBPmanual", "SBPauto"),
  TotalChol = c("TC")
)

PGS_factor_names <- list(
  BMI = c("BMI", "Extreme BMI", "Obesity (class I)", "Obesity (class II)", 
          "Obesity (class III)", "Overweight"),
  Diabetes = c("T2D"),
  Alcohol = c("DC2", "Alcohol Consumption"),
  Depression = c("MDD"),
  Diet = c("DC2"),
  Education = c("EA (Lee, 2014)", "EA (Okbay, 2022)"),
  HDL = c("HDL"),
  HeartDisease = c("CAD"),
  Physical = c("LST", "MVPA"),
  SysBP = c("SBP (manual reading)", "SBP (automated reading)"),
  TotalChol = c("TC")
)

factors <- c("BMI", "Diabetes", "Alcohol","Depression", "Diet", "Education",
             "HDL", "Physical","HeartDisease", "SysBP", "TotalChol")

factorNames <- c("BMI", "T2D", "L-M Alcohol","MDD", "Diet", "Education",
                 "HDL", "Physical Act.","HD", "SBP", "TC")

models <- names(df_list)
modelNames <- c("Lasso", "Lasso-sparse", "Ridge", 
                 "Bolt", "BayesR", "BayesR-shrink")

AUC_df <- NULL
for (i in 1:length(factors)){
  for (m in 1:length(models)){
    for (f in 1:length(PGS_factors[[factors[i]]])){
      samples1 <- rownames(Y_all[samples,])[!is.na(Y_all[samples,factors[i]])]
      roc_list <- roc(response = factor(Y_all[samples1, factors[i]], levels = c("No", "Yes")), 
                      predictor = df_list[[models[m]]][samples1,PGS_factors[[factors[i]]][f]])
      
      temp <- data.frame(AUC = auc(roc_list),
                           PGS = PGS_factor_names[[factors[i]]][f],
                           Model = modelNames[m],
                           Factor = factorNames[i])
      AUC_df <- rbind.data.frame(AUC_df, temp)
    }
   
  }
}

AUC_df$PGS <- factor(AUC_df$PGS, levels = unique(AUC_df$PGS))

colors <- c(RColorBrewer::brewer.pal(n = 8,"Dark2"),
            RColorBrewer::brewer.pal(n = 6,"Set2"))
p <- ggplot(AUC_df) +
  geom_bar(aes(x = PGS, y = AUC, fill = Factor, alpha = Model),
           stat = "identity", position = position_dodge()) +
  facet_grid(cols = vars(Factor), scales = "free", space = "free") +
  coord_cartesian(ylim = c(0.5,1)) +
  scale_alpha_manual(values = c(0.5,0.6,0.7,0.8,0.9,1)) +
  scale_fill_manual(values = colors) +
  guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")

ggsave(p, file = "PGS/PGS_prediction.png", width = 10, height = 6)
