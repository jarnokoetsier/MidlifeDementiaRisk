# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(ranger)
library(caret)
library(foreach)
library(doParallel)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(pROC)


factorNames <- c("Education", "Systolic BP", "BMI", "Total Chol.", "Physical Act.",
                 "Diet", "Smoking", "L-M Alcohol", "Depression", "Diabetes",
                 "HDL Chol.", "Heart Disease", "APOE", "Sex", "Age < 47", "Age > 53")

plotDF <- NULL
methodNames <- c("ElasticNet", "Random Forest")
for (m in 1:length(methodNames)){
  
  load(paste0("PerFactor/",methodNames[m],"/finalOutput_test.RData"))
  load(paste0("PerFactor/",methodNames[m],"/finalOutput_CV.RData"))
  
  AUC_all <- rep(NA,length(finalOutput_test))
  for (f in 1:length(finalOutput_test)){
    AUC_all[f] <- finalOutput_test[[f]]$AUC
  }
  temp_test <- data.frame(AUC = AUC_all,
                          Factor = factorNames,
                          Method = rep(methodNames[m], length(AUC_all)),
                          Set = rep("Test", length(AUC_all)))
  
  
  AUC_all <- rep(NA,length(finalOutput))
  for (f in 1:length(finalOutput)){
    AUC_all[f] <- finalOutput[[f]]$AUC
  }
  temp_CV <- data.frame(AUC = AUC_all,
                        Factor = factorNames,
                        Method = rep(methodNames[m], length(AUC_all)),
                        Set = rep("Cross-validation", length(AUC_all)))
  
  temp_all <- rbind.data.frame(temp_test, temp_CV)
  plotDF <- rbind.data.frame(plotDF, temp_all)
}


p <- ggplot(plotDF) +
  geom_bar(aes(x = Factor, y = AUC, fill = Factor, alpha = Set), 
           stat = "identity", position=position_dodge()) +
  facet_grid(rows = vars(Method)) +
  coord_cartesian(ylim = c(0.5,1)) +
  theme_bw()+
  guides(fill = "none") +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(n = 8, name = "Set2"),
                                RColorBrewer::brewer.pal(n = 8, name = "Dark2"))) +
  scale_alpha_manual(values = c(0.5,1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.title = element_blank())

ggsave(p, file = "PerFactor/AUC_factors.png", width = 10, height = 6)



n_cases <- rep(NA, length(finalOutput_test))
n_control <- rep(NA, length(finalOutput_test))
for (i in 1:length(finalOutput_test)){
  n_cases[i] <- sum(finalOutput_test[[i]]$ObsPred_test$obs == "Yes")
  n_control[i] <- sum(finalOutput_test[[i]]$ObsPred_test$obs == "No")
}

plotDF_n <- data.frame(Number = c(n_cases,n_control),
                       CaseControl = c(rep("Case", length(n_cases)),
                                    rep("Control", length(n_control))),
                       Factor = factorNames)

p <- ggplot(plotDF_n)+
  geom_bar(aes(x = Factor, y = Number, fill = CaseControl), 
           stat = "identity", color = "black") +
  coord_flip() +
  ylab("Number of samples") +
  xlab(NULL) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

ggsave(p, file = "PerFactor/CaseControlDistribution.png", width = 8, height = 6)