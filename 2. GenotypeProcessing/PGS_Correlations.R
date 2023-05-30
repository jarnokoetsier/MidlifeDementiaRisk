# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(patchwork)
library(data.table)
library(corrr)
library(ggdendroplot)

# Load PGSs
load("E:/Thesis/EXTEND/df_list.RData")
PGSresult <- df_list$bayesr

# Change names of PGSs
PGS_names <- c("AD (w/o APOE)", "AD (stage I)", "AD (stage I-II)", "Adiponectin", "Alchol consumption",
               "ASD", "BD", "BL", "BMI", "BW", "BW (fetal genetics)", "BW (maternal genetics)",
               "CAD", "Childhood obsesity", "Chronotype", "CKD (European)", "CKD (trans-ethnic)",
               "DC2", "EA (Lee, 2014)", "EA (Okbay, 2022)","Extreme BMI", "Extreme height",
               "Extreme HWR", "AD (family history)", "HbA1c", "HC", "HDL", "Infant Head Circumference",
               "LDL", "LST","MDD", "MVPA","Obesity (class I)", "Obsesity (Class II)", "Obesity (Class III)",
               "Overweigth", "PD", "Pubertal Growth", "RA", "SBP (automated reading)", "SBP (manual reading)",
               "SHR", "Sleep duration", "T2D", "TC",
               "TG", "WC", "WHR")

colnames(PGSresult) <- PGS_names

#*****************************************************************************#
#   PGS-PGS correlations
#*****************************************************************************#

# Calculate correlations
corrDF <- as.data.frame(correlate(PGSresult, diagonal = 1, method = "spearman"))
rownames(corrDF) <- corrDF$term
corrDF <- corrDF[,-1]

# Format data for plotting
plotCor <- gather(corrDF)
plotCor$key1 <- rep(rownames(corrDF), ncol(corrDF))

# Perform clustering to get sample order
model <- hclust(as.dist(1-abs(corrDF)), "ward.D2")
order <- model$labels[model$order]
plotCor$key <- factor(plotCor$key, levels = order)
plotCor$key1 <- factor(plotCor$key1, levels = order)

# Main correlation plot
main <- ggplot(plotCor) +
  geom_tile(aes(x = key, y = key1, fill = value),
            width = 0.95, height = 0.95) +
  scale_fill_gradient2(low = "#000072", mid = "#FFFBF5", high = "red", midpoint = 0,
                       limits = c(-1,1)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = paste0(str_to_title("spearman"),"\nCorrelation"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
        axis.text.y = element_text(size = 11))

# Dendrogram
dendroPlot <- ggplot() +
  geom_tile(data = plotCor, aes(x = as.numeric(key), y = 1), fill = "white", alpha = 0) +
  geom_dendro(model,xlim = c(1,length(unique(plotCor$key)))) +
  theme_void() +
  theme(legend.position = "none") 

# Combine plots
p <- dendroPlot + main +
  plot_layout(nrow = 2, ncol = 1,
              heights = c(1,6))

# Save plots
setwd("E:/Thesis/EXTEND/Phenotypes")
ggsave(p, file = "PGS_correlations.png", width = 11, height = 10)



#*****************************************************************************#
#   PGS-Phenotype correlations
#*****************************************************************************#

# samples
famFile <- fread("E:/Thesis/EXTEND/Genotypes/ChrBPData/Output_all/FINAL/EXTEND_PostImpute_FINAL_bp_dup.fam")

# Load meta data
setwd("E:/Thesis/EXTEND/Phenotypes")
load("metaData_ageFil.Rdata")

# High Education
Education <- ifelse(dat$None.of.the.above == 1,1,0)

# High Systolic Blood Pressure
SysBP <- ifelse(dat$MeanSysBP <= 140,0,1)

# High BMI
BMI <-ifelse(dat$BMI <= 30,0,1)

# High total Cholesterol
TotalChol <- ifelse(as.numeric(dat$Chol_unloged) <= 6.5,0,1)

# Low physical activity
Physical <- ifelse(dat$Exercise.increased.pulse.more.than.2halfhrsawk == 1,1,0)

# Healthy Diet
Diet <- ifelse(as.numeric(dat$Fruit) > 3 | as.numeric(dat$Vegtables) > 3, 1, 0)

# Smoking
Smoking <- ifelse(dat$Do.you.currently.smoke==1,1,0)

# Low alcohol intake
Alcohol <- rep(NA,nrow(dat))
Alcohol[dat$How.often.do.you.drink.alcohol==0] <- 1
Alcohol[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==1]<- 1
Alcohol[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==1]<-1
Alcohol[dat$alcoholic.drinks.per.day==3 &    dat$How.often.do.you.drink.alcohol==1]<-1
Alcohol[dat$alcoholic.drinks.per.day==4 &    dat$How.often.do.you.drink.alcohol==1]<-1
Alcohol[dat$alcoholic.drinks.per.day==5 &    dat$How.often.do.you.drink.alcohol==1]<-1

Alcohol[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==2]<-1
Alcohol[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==2]<-1

Alcohol[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==3]<-1
Alcohol[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==3]<-1

# High alcohol intake
Alcohol[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==1]<-0
Alcohol[dat$alcoholic.drinks.per.day==3  &  dat$How.often.do.you.drink.alcohol==2]<-0
Alcohol[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==2]<-0
Alcohol[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==2]<-0
Alcohol[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==3]<-0
Alcohol[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==3]<-0
Alcohol[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==3]<-0
Alcohol[dat$alcoholic.drinks.per.day==2 &   dat$How.often.do.you.drink.alcohol==4]<-0
Alcohol[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==4]<-0
Alcohol[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==4]<-0
Alcohol[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==4]<-0
Alcohol[dat$alcoholic.drinks.per.day=="." &    dat$How.often.do.you.drink.alcohol==4]<-0

# Depression
Depression <- ifelse(dat$Depression==1,1,0)

# Diabetes
Diabetes <-ifelse(dat$T2.Diabetes==1,1,0)

# High HDL
HDL <- ifelse(dat$HDL_unloged > 2.2, 1,0)

# Heart Disease
HeartDisease <- ifelse(dat$Heart.Disease==1,1,0)

# Kidney Disease
KidneyDisease <- ifelse(dat$Kidney..Disease==1,1,0)

# Age
Age<- rep(NA,nrow(dat))
Age[dat$Age < 47] <- 0
Age[dat$Age >= 47 & dat$Age <= 53] <- 3
Age[dat$Age > 53] <- 4

# Sex
Sex <- ifelse(dat$Sex == 1,1,0)

# APOE
APOEstatus <- read.csv("APOE_status.csv")
rownames(APOEstatus) <- APOEstatus$sampleID
APOEstatus <- APOEstatus[dat$ID,]
APOE <- ifelse(APOEstatus$e4 == 0,0,1)

# Combine factors
Y_all <- data.frame(Education,
                    SysBP,
                    BMI,
                    TotalChol,
                    Physical,
                    Diet,
                    Smoking,
                    Alcohol,
                    Depression,
                    Diabetes,
                    HDL,
                    HeartDisease,
                    KidneyDisease,
                    Age,
                    Sex,
                    APOE
)

rownames(Y_all) <- dat$ID

# Select samples
Y_PGS <- Y_all[famFile$V2,]

# Change names of PGSs for plotting
PGS_names <- c("AD (w/o APOE)", "AD (stage I)", "AD (stage I-II)", "Adiponectin", "Alchol consumption",
               "ASD", "BD", "BL", "BMI", "BW", "BW (fetal genetics)", "BW (maternal genetics)",
               "CAD", "Childhood obsesity", "Chronotype", "CKD (European)", "CKD (trans-ethnic)",
               "DC2", "EA (Lee, 2014)", "EA (Okbay, 2022)","Extreme BMI", "Extreme height",
               "Extreme HWR", "AD (family history)", "HbA1c", "HC", "HDL", "Infant Head Circumference",
               "LDL", "LST","MDD", "MVPA","Obesity (class I)", "Obesity (Class II)", "Obesity (Class III)",
               "Overweigth", "PD", "Pubertal Growth", "RA", "SBP (automated reading)", "SBP (manual reading)",
               "SHR", "Sleep duration", "T2D", "TC",
               "TG", "WC", "WHR")

# Calculate correlations
plotCor_all <- NULL
plotSig_all <- NULL
for (m in 1:length(df_list)){

  corMatrix <- matrix(NA, nrow = ncol(df_list[[m]]), ncol = ncol(Y_PGS)) # correlation coefficient
  sigMatrix <- matrix(NA, nrow = ncol(df_list[[m]]), ncol = ncol(Y_PGS)) # correlation significance
  for (i in 1:ncol(Y_PGS)){
    sigMatrix[,i] <- apply(df_list[[m]],2,function(x){cor.test(x,Y_PGS[,i], method = "spearman",
                                                               use = "pairwise.complete.obs")$p.value})
    corMatrix[,i] <- apply(df_list[[m]],2,function(x){cor.test(x,Y_PGS[,i], method = "spearman",
                                                               use = "pairwise.complete.obs")$estimate})
  }
  colnames(corMatrix) <- colnames(Y_PGS)
  colnames(corMatrix) <- c("Education", "Syst. BP", "BMI", "Total Chol." , "Physical Inact.", "Healthy Diet", "Smoking",
                           "L-M Alcohol", "Depression", "Type II Diabetes", "HDL Chol.", "Heart Disease", "Kidney Disease",
                           "Age", "Sex", "APOE \u03b54")
  rownames(corMatrix) <- PGS_names
  
  # Format data for plotting
  plotCor <- gather(as.data.frame(corMatrix))
  plotCor$key1 <- rep(rownames(corMatrix), ncol(corMatrix))
  plotCor$Method <- rep(names(df_list)[m], nrow(plotCor))
  plotCor_all <- rbind.data.frame(plotCor_all, plotCor)
  
  # Change column names
  colnames(sigMatrix) <- colnames(Y_PGS)
  colnames(sigMatrix) <- c("Education", "Syst. BP", "BMI", "Total Chol." , "Physical Inact.", "Healthy Diet", "Smoking",
                           "L-M Alcohol", "Depression", "Type II Diabetes", "HDL Chol.", "Heart Disease", "Kidney Disease",
                           "Age", "Sex", "APOE \u03b54")
  rownames(sigMatrix) <- colnames(df_list[[m]])
  
  # Format data for plotting
  plotSig <- gather(as.data.frame(sigMatrix))
  plotSig$key1 <- rep(rownames(sigMatrix), ncol(sigMatrix))
  plotSig$Method <- rep(names(df_list)[m], nrow(plotSig))
  
  plotSig_all <- rbind.data.frame(plotSig_all, plotSig)
}

plotCor_all$pvalue <- plotSig_all$value
plotCor_all$FDR <- p.adjust(plotCor_all$pvalue, method = "fdr")
plotCor_all$Sig <- ifelse(plotCor_all$FDR < 0.05, "Yes", "No")


# Main correlation plot
main <- ggplot(plotCor_all[plotCor_all$Method == "bayesr-shrink",]) +
  geom_tile(aes(x = key, y = key1, fill = value, color = Sig), 
            linewidth = 0.7, width = 0.8, height = 0.8) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.5,0.5)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = paste0(str_to_title("spearman"),"\nCorrelation"))+
  scale_color_manual(values = c("white","black")) +
  guides(color = "none") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  #theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))


# Side plot: to which score belongs the factor
whichScore <- data.frame(
  Factor <- c(c("Age", "Sex", "BMI", "Education", "Total Chol.", 
              "Physical Inact.", "Syst. BP"),
              c("Age", "Sex", "APOE \u03b54", "BMI", "Education", "Total Chol.", 
                "Physical Inact.", "Syst. BP"),
              c("Syst. BP", "BMI",  "Physical Inact.", "Healthy Diet", "Smoking",
                "L-M Alcohol", "Depression", "Type II Diabetes", "HDL Chol.", 
                "Heart Disease", "Kidney Disease")),
  Score <- c(rep("CAIDE1",7), rep("CAIDE2", 8), rep("LIBRA",11))
)

bottom <- ggplot(whichScore) +
  geom_tile(aes(x = Factor, y = Score, fill = Score), 
            color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_manual(values = c("#DC3535","#F49D1A","#B01E68")) +
  theme_bw() +
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#DCDCDC"),
        legend.position = "none")

# Combine plots
p <- main + bottom + 
plot_layout(nrow = 2, ncol = 1,
            heights = c(10,1))
# Save plot
setwd("E:/Thesis/EXTEND/Phenotypes")
ggsave(p, file = "PRS_Score_Cor.png", width = 8, height = 10)


# Horizontal Version

# Main correlation plot
main <- ggplot(plotCor_all[plotCor_all$Method == "bayesr-shrink",]) +
  geom_tile(aes(x = key1, y = key, fill = value, color = Sig), 
            linewidth = 0.7, width = 0.8, height = 0.9) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.5,0.5)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = paste0(str_to_title("spearman"),"\nCorrelation"))+
  scale_color_manual(values = c("white","black")) +
  guides(color = "none") +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Side plot: to which score belongs the factor
whichScore <- data.frame(
  Factor <- c(c("Age", "Sex", "BMI", "Education", "Total Chol.", 
                "Physical Inact.", "Syst. BP"),
              c("Age", "Sex", "APOE \u03b54", "BMI", "Education", "Total Chol.", 
                "Physical Inact.", "Syst. BP"),
              c("Syst. BP", "BMI",  "Physical Inact.", "Healthy Diet", "Smoking",
                "L-M Alcohol", "Depression", "Type II Diabetes", "HDL Chol.", 
                "Heart Disease", "Kidney Disease")),
  Score <- c(rep("CAIDE1",7), rep("CAIDE2", 8), rep("LIBRA",11))
)

bottom <- ggplot(whichScore) +
  geom_tile(aes(y = Factor, x = Score, fill = Score), 
            color = "black") +
  xlab(NULL) +
  ylab(NULL) +
  scale_fill_manual(values = c("#DC3535","#F49D1A","#B01E68")) +
  theme_bw() +
  theme(axis.text.y = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#DCDCDC"),
        legend.position = "none")

# Combine plots
p <- bottom + main+
  plot_layout(nrow = 1, ncol = 2,
              widths = c(1,8))

# Save plots
setwd("E:/Thesis/EXTEND/Phenotypes")
ggsave(p, file = "PRS_Score_Cor_horizontal.png", width = 10, height = 6)
