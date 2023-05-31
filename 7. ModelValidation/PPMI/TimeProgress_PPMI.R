# Load packages
library(tidyverse)
library(caret)
library(patchwork)
library(ranger)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
cogcat <- read.csv("PPMI/CogCatWide_Filtered.csv")
cogcat$PATNO <- as.character(cogcat$PATNO)
load("PPMI/predictedScore_factors_PPMI.RData")
load("PPMI/metaData_ppmi.RData")

# Combine scores and meta data
metaData_all <- metaData_all[(metaData_all$age >= 40) & (metaData_all$age <= 75),]
samples <- intersect(metaData_all$Basename, rownames(predictedScore_factors))
rownames(metaData_all) <- metaData_all$Basename
metaData_fil <- metaData_all[samples,]
predictedScore_factors_fil <- predictedScore_factors[samples,]
length(unique(metaData_fil$ID))

###############################################################################

# Proportion of cognitive impairment at each time point

###############################################################################

# Get cognitive status at each time point
test <- inner_join(metaData_fil, cogcat, by = c("PATNO" = "PATNO"))
table(test$Class, test$CogDecon)
predictedScore_factors_fil <- predictedScore_factors_fil[test$Basename,]

# Make low and high risk category
load("~/PPMI/Fit_EMIF_MCI_RF.RData")
pred_RF <- predict(fit, predictedScore_factors_fil, type = "prob")
quantile(pred_RF$MCI,0.5)

test1 <- test[pred_RF$MCI < quantile(pred_RF$MCI,0.5),]
test2 <- test[pred_RF$MCI >= quantile(pred_RF$MCI,0.5),] 

# Format data for plotting
n_normal_low <- rep(NA, 9)
n_total_low <- rep(NA, 9)
ci_up_low <- rep(NA, 9)
ci_down_low <- rep(NA, 9)
ratio_low <- rep(NA, 9)
n_normal_high <- rep(NA, 9)
n_total_high <- rep(NA, 9)
ci_up_high <- rep(NA, 9)
ci_down_high <- rep(NA, 9)
ratio_high <- rep(NA, 9)
for(j in 1:9){
  i = j-1
  
  n_normal_low[j] <- sum(test1[!is.na(test1[,paste0("X",i)]), paste0("X",i)] == "Normal" | 
                    test1[!is.na(test1[,paste0("X",i)]), paste0("X",i)] == "Cognitive Complaint")
  n_total_low[j] <- sum(!is.na(test1[,paste0("X",i)]))
  ratio_low[j] <- n_normal_low[j]/n_total_low[j]
  ci_up_low[j] <- binom.test(n_normal_low[j], n_total_low[j])$conf.int[[2]]
  ci_down_low[j] <- binom.test(n_normal_low[j], n_total_low[j])$conf.int[[1]]
  
  
  n_normal_high[j] <- sum(test2[!is.na(test2[,paste0("X",i)]), paste0("X",i)] == "Normal" | 
                           test2[!is.na(test2[,paste0("X",i)]), paste0("X",i)] == "Cognitive Complaint")
  n_total_high[j] <- sum(!is.na(test2[,paste0("X",i)]))
  ratio_high[j] <- n_normal_high[j]/n_total_high[j]
  ci_up_high[j] <- binom.test(n_normal_high[j], n_total_high[j])$conf.int[[2]]
  ci_down_high[j] <- binom.test(n_normal_high[j], n_total_low[j])$conf.int[[1]]
}

# Combine into data frame
plotDF <- data.frame(n_normal_low,
                     n_total_low,
                     ci_up_low,
                     ci_down_low,
                     ratio_low,
                     n_normal_high,
                     n_total_high,
                     ci_up_high,
                     ci_down_high,
                     ratio_high,
                     Time = 0:8)

# Test for significant difference at each time point
test <- plotDF[plotDF$Time == 0,]
m <- matrix(c(test$n_normal_low, test$n_total_low - test$n_normal_low,
       test$n_normal_high, test$n_total_high - test$n_normal_high),
       nrow = 2, ncol = 2)
chisq.test(m)


# Make plot
main <- ggplot(plotDF) +
  geom_point(aes(x = Time, y = ratio_low, 
                 color = "Low Risk")) +
  geom_step(aes(x = Time, y = ratio_low, 
                color = "Low Risk")) +
  geom_point(aes(x = Time, y = ratio_high,
                 color = "High Risk")) +
  geom_step(aes(x = Time, y = ratio_high,
                color = "High Risk")) +
  geom_bar(aes(x = Time, y = n_total_low/80), 
           stat = "identity", position = position_dodge(), alpha = 0) +
  ylim(c(0.3,1)) +
  ylab("Probability of\nnormal cognition") +
  xlab("Time (years)") +
  scale_color_manual(values = c("#CB181D","#2171B5")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

top <- ggplot(plotDF) +
  geom_bar(aes(x = Time, y = n_total_low, fill = "Low Risk"), 
           stat = "identity", position = position_dodge(), color = "black") +
  geom_bar(aes(x = Time, y = n_total_high, fill = "High Risk"),
           stat = "identity", position = position_dodge(), width = 0.7, color = "black") +
  xlab(NULL) +
  ylab("# Samples") +
  scale_fill_manual(values = c("#CB181D","#2171B5")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.position = "none")

p <- top + main +
  plot_layout(ncol = 1, nrow = 2,
              heights = c(1, 5))

# Save plot
ggsave(p, file = "TimeAnalysis_PPMI.png", width = 6, height = 6)



###############################################################################

# Kaplan-Meier

###############################################################################


# Load packages
library(tidyverse)
library(caret)
library(patchwork)
library(ranger)
library(survival)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(tidyverse)
library(tidyverse)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
cogcat <- read.csv("PPMI/CogCatWide_Filtered.csv")
cogcat$PATNO <- as.character(cogcat$PATNO)
load("PPMI/predictedScore_factors_PPMI.RData")
load("PPMI/metaData_ppmi.RData")

# Remove reverters
reverters <- cogcat
reverters[is.na(reverters)] <- ""

rev_cat <- ""
for (i in 3:11){
  rev_cat <- paste0(rev_cat, reverters[,i]) 
}

reverters <- c(cogcat$PATNO[str_detect(rev_cat, "DementiaNormal")],
               cogcat$PATNO[str_detect(rev_cat, "MCINormal")],
               cogcat$PATNO[str_detect(rev_cat, "DementiaCognitive Complaint")],
               cogcat$PATNO[str_detect(rev_cat, "MCICognitive Complaint")])

cogcat <- cogcat[!(cogcat$PATNO %in% reverters),]

# midlife samples
metaData_all <- metaData_all[(metaData_all$age >= 40) & (metaData_all$age <= 75),]
samples <- intersect(metaData_all$Basename, rownames(predictedScore_factors))
rownames(metaData_all) <- metaData_all$Basename
metaData_fil <- metaData_all[samples,]
predictedScore_factors_fil <- predictedScore_factors[samples,]

# combine meta data with cogcat
test <- inner_join(metaData_fil, cogcat, by = c("PATNO" = "PATNO"))
table(test$Class, test$CogDecon)
predictedScore_factors_fil <- predictedScore_factors_fil[test$Basename,]

# make predictions
load("~/PPMI/Fit_EMIF_MCI_RF.RData")
pred_RF <- predict(fit, predictedScore_factors_fil, type = "prob")

predictDF <- data.frame(PATNO = test$PATNO, 
                        pred = pred_RF$MCI)

# Split into low and high risk class
predictDF$predClass <- ifelse(predictDF$pred < quantile(pred_RF$MCI,0.5),"Low MCI Risk (RF)","High MCI Risk (RF)")


# Prepare data for analysis
testDF <- gather(as.data.frame(cogcat[,3:11]))
testDF$PATNO <- rep(cogcat$PATNO, 9)
testDF <- testDF[!is.na(testDF$value),]
testDF <- testDF[(testDF$value == "Dementia") | (testDF$value == "MCI"),]
testDF$Time <- as.numeric(str_remove(testDF$key, "X"))
testDF$Status <- ifelse((testDF$value == "Dementia") | (testDF$value == "MCI"), "Cognitive Impaired", "Normal")

for (i in unique(testDF$PATNO)){
  testDF[testDF$PATNO == i, "Time"] <- min(testDF[testDF$PATNO == i, "Time"])
}

testDF <- testDF[!duplicated(testDF[,3:5]),]
normal <- data.frame(key = "X9",
                     value = "Normal",
                     PATNO = setdiff(cogcat$PATNO, unique(testDF$PATNO)),
                     Time = 9,
                     Status = "Normal")
testDF <- rbind.data.frame(testDF, normal)
length(unique(cogcat$PATNO))


kaplanDF <- testDF[,c("PATNO", "Time", "Status")]
# 1: censored
# 2: disease
kaplanDF$Test <- ifelse(kaplanDF$Status == "Normal",1,2)
kaplanDF <- inner_join(kaplanDF, predictDF, by = c("PATNO" = "PATNO"))


# Perform time analysis
test <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF)

# Prepare data for plotting
plotDF <- data.frame(
  Time = c(test$time),
  n.risk = c(test$n.risk),
  n.event = c(test$n.event),
  n.censor = c(test$n.censor),
  surv = c(test$surv),
  Class = factor(c(rep("High MCI Risk (RF)",10), rep("Low MCI Risk (RF)",6)),
                 levels = c("Low MCI Risk (RF)", "High MCI Risk (RF)")))


# Make Kaplan-meier curve
survdiff(Surv(Time, Test) ~ predClass, data = kaplanDF)
kaplanDF$predClass <- factor(kaplanDF$predClass, levels = c("Low MCI Risk (RF)", "High MCI Risk (RF)"))

p <- survfit2(Surv(Time, Test) ~ predClass, data = kaplanDF) %>% 
  ggsurvfit(size = 1.5) +
  add_confidence_interval() +
  scale_color_manual(values = c("#FD8D3C","#D94801")) +
  scale_fill_manual(values = c("#FD8D3C","#D94801")) +
  theme_classic() +
  ylab("Probability of\nnormal cognition") +
  xlab("Time (years)") +
  scale_x_continuous(breaks = c(0,2,4,6,8)) +
  #ggtitle("Cognitive impairment (PPMI)") +
  #xlim(c(0,8)) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))

# Save plot
ggsave(p, file = "KaplanMeier_PPMI_MRSonly.png", width = 8, height = 4.5)
