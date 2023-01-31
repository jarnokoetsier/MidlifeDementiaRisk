library(tidyverse)
library(patchwork)
setwd("E:/Thesis/EXTEND/Genotypes")
load("df_Result_PGS.RData")

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


Y_PGS <- Y_all[rownames(df_Result_PGS),]


corMatrix <- matrix(NA, nrow = ncol(df_Result_PGS), ncol = ncol(Y_PGS))
for (i in 1:ncol(Y_PGS)){
  corMatrix[,i] <- apply(df_Result_PGS,2,function(x){cor(x,Y_PGS[,i], method = "spearman",
                                                   use = "pairwise.complete.obs")})
}


colnames(corMatrix) <- colnames(Y_PGS)
colnames(corMatrix) <- c("Education", "Syst. BP", "BMI", "Total Chol." , "Physical Act.", "Healthy Diet", "Smoking",
  "L-M Alcohol", "Depression", "Type II Diabetes", "HDL Chol.", "Heart Disease", "Kidney Disease",
  "Age", "Sex", "APOE \u03b54")
rownames(corMatrix) <- colnames(df_Result_PGS)

# Format data for plotting
plotCor <- gather(as.data.frame(corMatrix))
plotCor$key1 <- rep(rownames(corMatrix), ncol(corMatrix))



################################################################################
sigMatrix <- matrix(NA, nrow = ncol(df_Result_PGS), ncol = ncol(Y_PGS))
for (i in 1:ncol(Y_PGS)){
  sigMatrix[,i] <- apply(df_Result_PGS,2,function(x){cor.test(x,Y_PGS[,i], method = "spearman",
                                                         use = "pairwise.complete.obs")$p.value})
}


colnames(sigMatrix) <- colnames(Y_PGS)
colnames(sigMatrix) <- c("Education", "Syst. BP", "BMI", "Total Chol." , "Physical Act.", "Healthy Diet", "Smoking",
                         "L-M Alcohol", "Depression", "Type II Diabetes", "HDL Chol.", "Heart Disease", "Kidney Disease",
                         "Age", "Sex", "APOE \u03b54")
rownames(sigMatrix) <- colnames(df_Result_PGS)

plotSig <- gather(as.data.frame(sigMatrix))
plotSig$key1 <- rep(rownames(sigMatrix), ncol(sigMatrix))
plotSig$FDR <- p.adjust(plotSig$value, method = "fdr")
plotSig_fil <- plotSig[plotSig$FDR < 0.05,]

plotCor$Sig <- rep("No", nrow(plotCor))
for (i in 1:nrow(plotSig_fil)){
  plotCor$Sig[(plotCor$key == plotSig_fil$key[i]) &(plotCor$key1 == plotSig_fil$key1[i])] <- "Yes"
}

################################################################################


# Highest correlation for each factor
maxCor <- data.frame(
  Factor = names(apply(corMatrix,2,function(x){which.max(abs(x))})),
  MaxPGS = rownames(corMatrix)[apply(corMatrix,2,function(x){which.max(abs(x))})])

plotCor <- inner_join(plotCor, maxCor, by = c("key" =  "Factor"))
plotCor$MaxPGS[plotCor$MaxPGS != plotCor$key1] <- "No"
plotCor$MaxPGS[plotCor$MaxPGS == plotCor$key1] <- "Yes"

# Main correlation plot
main <- ggplot(plotCor) +
  geom_tile(aes(x = key, y = key1, fill = value, color = Sig), 
            linewidth = 0.7, width = 0.8, height = 0.9) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.5,0.5)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = paste0(str_to_title("spearman"),"\nCorrelation"))+
  scale_color_manual(values = c("white","black")) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())



whichScore <- data.frame(
  Factor <- c(c("Age", "Sex", "BMI", "Education", "Total Chol.", 
              "Physical Act.", "Syst. BP"),
              c("Age", "Sex", "APOE \u03b54", "BMI", "Education", "Total Chol.", 
                "Physical Act.", "Syst. BP"),
              c("Syst. BP", "BMI",  "Physical Act.", "Healthy Diet", "Smoking",
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

p <- main + bottom + 
plot_layout(nrow = 2, ncol = 1,
            heights = c(8,1))

setwd("E:/Thesis/EXTEND/Phenotypes")
ggsave(p, file = "PRS_Score_Cor.png", width = 8, height = 10)


# Horizontal

# Main correlation plot
main <- ggplot(plotCor) +
  geom_tile(aes(x = key1, y = key, fill = value, color = Sig), 
            linewidth = 0.7, width = 0.8, height = 0.9) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.5,0.5)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = paste0(str_to_title("spearman"),"\nCorrelation"))+
  scale_color_manual(values = c("white","black")) +
  guides(color = "none") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



whichScore <- data.frame(
  Factor <- c(c("Age", "Sex", "BMI", "Education", "Total Chol.", 
                "Physical Act.", "Syst. BP"),
              c("Age", "Sex", "APOE \u03b54", "BMI", "Education", "Total Chol.", 
                "Physical Act.", "Syst. BP"),
              c("Syst. BP", "BMI",  "Physical Act.", "Healthy Diet", "Smoking",
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

p <- bottom + main+
  plot_layout(nrow = 1, ncol = 2,
              widths = c(1,8))

setwd("E:/Thesis/EXTEND/Phenotypes")
ggsave(p, file = "PRS_Score_Cor_horizontal.png", width = 10, height = 6)
