
# Load packages
library(readxl)
library(tidyverse)

# Set working directory
setwd("E:/Thesis/EXTEND/Phenotypes")

# Load data
metadata<-read.csv("extend_pheno_all.csv")
meth_ID<-read.csv("methylation_phenotype_oct2020.csv")

# Select subset
meth_ID <- meth_ID[,c("Basename","IID","Chip","Position","Plate")]   
colnames(meth_ID) <- c("Basename","ID","Chip","Position","Plate")
colnames(metadata)[2] <-"ID"

# Merge data frames
dat <- merge(meth_ID, metadata, all.x=T)

# Remove samples with missing age
dat <- dat[!is.na(dat$Age),]

# Select samples with age between 75 and 40 only
dat <- dat[dat$Age <= 75 & dat$Age >= 40,]
rm(metadata, meth_ID)

dat_copy <- dat

###############################################################################

# Calculate CAIDE1 score

###############################################################################

# Information about the calculation of the CAIDE scores can be found here:
# https://academic.oup.com/biomedgerontology/article/76/8/1407/6273220?login=true


#*****************************************************************************#
# Age
#*****************************************************************************#
# Discretize age in 3 classes with score 0, 3, and 4.
dat$age_c <- rep(NA,nrow(dat))
dat$age_c[dat$Age < 47] <- 0
dat$age_c[dat$Age >= 47 & dat$Age <= 53] <- 3
dat$age_c[dat$Age > 53] <- 4
table(dat$age_c)

#*****************************************************************************#
# Sex
#*****************************************************************************#
# Male = 1 (1 in dataset)
# Female = 0 (2 in dataset)

dat$Sex_c <- ifelse(dat$Sex == 1,1,0) # Check!!!
table(dat$Sex_c)


#*****************************************************************************#
# Education
#*****************************************************************************#
# "College.or.Uni.degree"                          "A.level.AS.level.or.equiv"                                    
# "O.level.GCSEs.or.equiv"                         "CSEs.or.equiv"                                                
# "NVQ.HND.HNC.or.equiv"                            "Other.professional.quals"                                     
# "None.of.the.above
dat$Edu_c <- ifelse(dat$None.of.the.above == 1,2,0)
table(dat$Edu_c)


#*****************************************************************************#
# Systolic blood pressure
#*****************************************************************************#
dat$Syst_c <- ifelse(dat$MeanSysBP <= 140,0,2)
table(dat$Syst_c)

#*****************************************************************************#
# BMI
#*****************************************************************************#
dat$BMI_c<-ifelse(dat$BMI <= 30,0,2)
table(dat$BMI_c)

#*****************************************************************************#
# Serum total cholesterol level
#*****************************************************************************#
dat <- dat[-which(dat$Chol_unloged=="."),]
dat$Chol_unloged <- as.numeric(dat$Chol_unloged)
dat$Chol_c <- ifelse(dat$Chol_unloged <= 6.5,0,2)
table(dat$Chol_c)

#*****************************************************************************#
# Physical activity
#*****************************************************************************#
dat$PHYSICAL_c <- ifelse(dat$Exercise.increased.pulse.more.than.2halfhrsawk == 1,0,1)
table(dat$PHYSICAL_c)


#*****************************************************************************#
# Calculate scores
#*****************************************************************************#
CAIDE <- dat[,c(1:8,11,12,151:156)]
CAIDE$X <- NULL

CAIDE$CAIDE <- rowSums(CAIDE[,10:15])
hist(CAIDE$CAIDE,10)

# Save data
CAIDE <- unique(CAIDE)
save(CAIDE, file="CAIDE.Rdata")

###############################################################################

# Calculate CAIDE2 score

###############################################################################

#APOE: The two SNPs (rs429358, rs7412) that define the epsilon 2, 3, and 4 alleles



###############################################################################

# Calculate LIBRA score

###############################################################################

dat <- dat_copy

#*****************************************************************************#
# Adherence to a Mediterranean diet (-1.7)
#*****************************************************************************#

# Reported amount of fruits and vegetables consumed by the participant the previous day. 
# A healthy diet was defined as consuming five or more portions of fruits and vegetables on a daily basis [8].
# [8]	National Health Service (2009) 5 A Day. NHS, London. (CHECK FOR UPDATE)
# Variables in EXTEND: 
# dat$Fruit     ==  Fruit(portions per day)	        0    1-2     3-4     5-6     >6
# dat$Vegetables ==  Vegetables(portions per day)	  0    1-2     3-4     5-6     >6

dat <- dat[-which(dat$Fruit=="."|dat$Vegtables=="."),]
dat$Fruit <- as.numeric(dat$Fruit)
dat$Vegtables < -as.numeric(dat$Vegtables)

dat$MEDITERANIAN <- ifelse(dat$Fruit > 3 | dat$Vegtables > 3, 1,0)
table(dat$MEDITERANIAN)

# Add weight to column
dat$MEDITERANIAN<-ifelse(dat$Fruit > 3 | dat$Vegtables > 3, -1.7,0)


#*****************************************************************************#
# Physical inactivity (+1.1)
#*****************************************************************************#

# Variables in EXTEND: Exercise causing  increased pulse more than 2.5 hours per week
dat <- dat[-which(dat$Exercise.increased.pulse.more.than.2halfhrsawk=="."),]
dat$PHYSICAL_INACTIVITY<-ifelse(dat$Exercise.increased.pulse.more.than.2halfhrsawk == 0, 1, 0)
table(dat$PHYSICAL_INACTIVITY)

# Add weight to column
dat$PHYSICAL_INACTIVITY<-ifelse(dat$Exercise.increased.pulse.more.than.2halfhrsawk == 0, 1.1,0)


#*****************************************************************************#
# Smoking (+1.5)
#*****************************************************************************#

# Variables in EXTEND: Do you currently smoke
dat <- dat[!is.na(dat$Do.you.currently.smoke),]
dat$SMOKING <- ifelse(dat$Do.you.currently.smoke==1,1,0)
table(dat$SMOKING)

# Add weight to column
dat$SMOKING<-ifelse(dat$Do.you.currently.smoke==1,1.5,0)


#*****************************************************************************#
# Low-to-moderate alcohol intake (-1)
#*****************************************************************************#

# UK: Low-to-moderate alcohol use was defined as 1-14 glasses per week according to recent UK alcohol guidelines [7].
#[7]	Department of Health (2016) Alcohol guidelines review - report from the guidelines development group to the UK Chief Medical Officers. Department of Health.
# Variables in EXTEND:
table(dat$alcoholic.drinks.per.day)
table(dat$How.often.do.you.drink.alcohol)
dat <- dat[-which(dat$How.often.do.you.drink.alcohol=="."),]
x <- dat[which(dat$alcoholic.drinks.per.day=="."),c("alcoholic.drinks.per.day","How.often.do.you.drink.alcohol")]

#alcoholic.drinks.per.day       ==1-2(1)     3-4(2)     5-6(3)    7-9(4)     >10(5)
#How.often.do.you.drink.alcohol ==Never(0)    1 a month(1)     2-4 mth(2)    2-3 wk(3)    =>4 wk(4)

# Low-moderate alcohol intake
dat$LtoMAlcohol <- NULL
dat[dat$How.often.do.you.drink.alcohol==0 ,"LtoMAlcohol"] <- 1
dat[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==3 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==4 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==5 &    dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-1

dat[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-1

dat[dat$alcoholic.drinks.per.day==1 &    dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-1
dat[dat$alcoholic.drinks.per.day==2 &    dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-1

# High alcohol intake
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==1 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==3  &  dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==2 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==3 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==2 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==3 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==4 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day==5 &   dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0
dat[dat$alcoholic.drinks.per.day=="." &    dat$How.often.do.you.drink.alcohol==4 ,"LtoMAlcohol"]<-0

# Remove NA rows
dat<-dat[!is.na(dat$LtoMAlcohol),]

# Check table
table(dat$LtoMAlcohol)

# Add weight to column
dat$LtoMAlcohol<-ifelse(dat$LtoMAlcohol == 1,-1,0)


#*****************************************************************************#
# Obesity (+1.6)
#*****************************************************************************#

dat$OBESITY <- ifelse(dat$BMI >=30, 1, 0)
table(dat$OBESITY)

# Add weight to column
dat$OBESITY<-ifelse(dat$BMI >=30,1.6,0)


#*****************************************************************************#
# Depression (+2.1)
#*****************************************************************************#

table(dat$Depression)
dat$DEPRESSION<-dat$Depression

# Add weight to column
dat$DEPRESSION <- ifelse(dat$DEPRESSION==1,2.1,0)


#*****************************************************************************#
# Type-2-Diabetes (+1.3)
#*****************************************************************************#

table(dat$T2.Diabetes)
dat$DIABETEII <- dat$T2.Diabetes

# Add weight to column
dat$DIABETEII<-ifelse(dat$DIABETEII==1,1.3,0)

#*****************************************************************************#
# Hypertension (+1.6)
#*****************************************************************************#
# Mean systolic blood pressure >= 140 mm Hg or mean diastolic blood pressure >= 90 mm Hg [3].

# No missing values
table(is.na(dat$MeanSysBP))
table(is.na(dat$MeanDiaBP))

dat$HYPERTENSTION <- ifelse (dat$MeanSysBP >=140 | dat$MeanDiaBP >= 90,1,0)
table(dat$HYPERTENSTION)

# Add weight to column
dat$HYPERTENSTION<-ifelse (dat$MeanSysBP >= 140 | dat$MeanDiaBP >= 90,1.6,0)

#*****************************************************************************#
# High cholesterol(+1.4)
#*****************************************************************************#

# Total cholesterol level of >= 5.0 mmol/L and low-density lipoprotein of >= 3.0 mmol/L, 
# following the guidelines of the National Health Service UK [2].

# Remove missing values
dat <- dat[!is.na(dat$HDL_unloged),]
dat$HDL_unloged <- as.numeric(dat$HDL_unloged)

#https://academic.oup.com/view-large/figure/342663278/cvab164f2.tif
dat$Highcholesterol<-ifelse(dat$HDL_unloged > 2.2, 1,0)
table(dat$Highcholesterol)

# Add weight to column
dat$Highcholesterol<-ifelse(dat$HDL_unloged > 2.2, 1.4,0)


#*****************************************************************************#
# Heart disease (+1.0)
#*****************************************************************************#

dat$Heartdisease <- dat$Heart.Disease

# Add weight to column
dat$Heartdisease<-ifelse(dat$Heartdisease==1,1,0)
table(dat$Heartdisease)

#*****************************************************************************#
# Chronic kidney disease (+1.1)
#*****************************************************************************#

# No missing values
table(is.na(dat$Kidney..Disease))
dat$kidneydisease<-dat$Kidney..Disease

# Add weight to column
dat$kidneydisease<-ifelse(dat$kidneydisease==1,1.1,0)
table(dat$kidneydisease)


#*****************************************************************************#
# Cognitive activity (-3.2)
#*****************************************************************************#

#NA

#*****************************************************************************#
# Calculate LIBRA scores
#*****************************************************************************#

EPILIBRA<-dat[,c(1:8,11,12,151:161)]
EPILIBRA$X<-NULL
EPILIBRA$LIBRA<-rowSums(EPILIBRA[,10:20])

# log scale
EPILIBRA$LIBRAlog10<- log10(EPILIBRA$LIBRA + 10)

# Save data
EPILIBRA <- unique(EPILIBRA)
save(EPILIBRA, file="EPILIBRA.Rdata")


###############################################################################

# Compare LIBRA and CAIDE

###############################################################################
library(ggpubr)
library(grid)

# Combine scores into data frame
scoreAll <- inner_join(CAIDE, EPILIBRA, by = c("ID" = "ID"))
scaleLIBRA <- (scoreAll$LIBRA - mean(scoreAll$LIBRA))/sd(scoreAll$LIBRA)
scaleCAIDE <- (scoreAll$CAIDE - mean(scoreAll$CAIDE))/sd(scoreAll$CAIDE)
scoreAll$RiskAll <- scaleLIBRA + scaleCAIDE

cor(scoreAll$LIBRA, scoreAll$CAIDE)

# Make scatter plot
scatter <- ggplot() +
  geom_point(data = scoreAll, aes(x = LIBRA, y = CAIDE, color = RiskAll),
             position=position_jitter(h=0.1,w=0.1), alpha = 0.5, size = 1.5) +
  #geom_smooth(data = scoreAll,method='lm', aes(x = LIBRA, y = CAIDE), formula = y ~x) +
  ylab("CAIDE1") +
  xlab("LIBRA") +
  labs(color = "") +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_continuous(breaks = c(-4,6), labels = c("Low risk", "High risk"),
                         type = "viridis")

# Make histogram of CAIDE scores
histogram_caide <- ggplot() +
  geom_histogram(data = scoreAll, aes(x = CAIDE), bins = 12) +
  coord_flip() +
  theme_void()

# Make histogram of LIBRA scores
histogram_libra <- ggplot() +
  geom_histogram(data = scoreAll, aes(x = LIBRA), bins = 12) +
  theme_void()

# Combine plots into single image
p <- ggarrange(histogram_libra,
          NULL,
          scatter,
          histogram_caide,
          nrow = 2,
          ncol = 2,
          widths = c(8,2),
          heights = c(2,8),
          align = "hv",
          common.legend = FALSE)

# Save plot
ggsave(p, file = "CAIDEvsLIBRA.png", height = 7, width = 8)


# Get legend
legendPlot <- ggplot() +
  geom_point(data = scoreAll, aes(x = LIBRA, y = CAIDE, color = RiskAll),
             position=position_jitter(h=0.1,w=0.1), alpha = 0.5, size = 1.5) +
  #geom_smooth(data = scoreAll,method='lm', aes(x = LIBRA, y = CAIDE), formula = y ~x) +
  labs(color = "") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_continuous(breaks = c(-4,6), labels = c("Low risk", "High risk"),
                         type = "viridis")
legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)


###############################################################################

# Save complete data

###############################################################################

dat <- dat_copy
save(dat, file = "metaData_ageFil.RData")


###############################################################################

# Overlap different scores

###############################################################################

# All samples (age filtered)
load("metaData_ageFil.RData")
all_samples <- dat$ID

# CAIDE1 samples
load("CAIDE.Rdata")
CAIDE1_samples <- CAIDE$ID

# CAIDE2 samples
load("E:/Thesis/EXTEND/Genotypes/Plots/gt_results.RData")
CAIDE2_samples <- names(gt_rs7412)[names(gt_rs7412) %in% CAIDE1_samples]

# LIBRA samples
load("EPILIBRA.Rdata")
LIBRA_samples <- EPILIBRA$ID

# Combine into single data frame
plotOverlap <- data.frame(ID = c(unique(all_samples),
                                 unique(CAIDE1_samples),
                                 unique(CAIDE2_samples),
                                 unique(LIBRA_samples)),
                          Source = c(rep("All",length(unique(all_samples))),
                                     rep("CAIDE1", length(unique(CAIDE1_samples))),
                                     rep("CAIDE2",length(unique(CAIDE2_samples))),
                                     rep("LIBRA",length(unique(LIBRA_samples))))
                          )

# determine order of samples and sources
orderID <- names(rev(sort(table(plotOverlap$ID[plotOverlap$Source != "LIBRA"]))))
orderID <- names(rev(sort(table(plotOverlap$ID))))
plotOverlap$ID <- factor(plotOverlap$ID, levels = orderID)
plotOverlap$Source <- factor(plotOverlap$Source,
                             c("All", "CAIDE1", "CAIDE2","LIBRA"))

# Make plot
p <- ggplot() +
  geom_tile(data = plotOverlap, aes(x = ID, y = Source, fill = Source)) +
  geom_text(aes(x = orderID[round(length(orderID)/2)],
                y =  c("All", "CAIDE1", "CAIDE2","LIBRA"),
                label = c(length(unique(all_samples)), length(unique(CAIDE1_samples)),
                          length(unique(CAIDE2_samples)), length(unique(LIBRA_samples)))),
            color = "white"
  ) +
  scale_fill_brewer(palette = "Set1") +
  xlab("Samples") +
  ylab("") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")

ggsave(p, file = "OverlapData.png", width = 8, height = 6)


###############################################################################

# Plot age distribution

###############################################################################
# Set working directory
setwd("E:/Thesis/EXTEND/Phenotypes")

# Load data
metadata<-read.csv("extend_pheno_all.csv")
meth_ID<-read.csv("methylation_phenotype_oct2020.csv")

# Select subset
meth_ID <- meth_ID[,c("Basename","IID","Chip","Position","Plate")]   
colnames(meth_ID) <- c("Basename","ID","Chip","Position","Plate")
colnames(metadata)[2] <-"ID"

# Merge data frames
dat <- merge(meth_ID, metadata, all.x=T)

# Remove samples with missing age
metaData <- dat[!is.na(dat$Age),]

# Prepare meta data object
metaData$AgeGroup <- rep("Mid-life",nrow(metaData))
metaData$AgeGroup[metaData$Age < 40] <- "Early-life"
metaData$AgeGroup[metaData$Age > 75] <- "Late-life"
metaData$AgeGroup <- factor(metaData$AgeGroup,
                            levels = c("Early-life",
                                       "Mid-life",
                                       "Late-life"))

# Make color palette
pal <- RColorBrewer::brewer.pal(n = 4, name = "Reds")[2:4]
names(pal) <- c("Early-life", "Mid-life", "Late-life")

x_min <- 40
x_max <- 75
# Make top of plot
top <- ggplot() +
  geom_rect(aes(xmin = min(metaData$Age)-5, xmax = x_min, ymin = 0, ymax = 1), fill = pal[1], alpha = 0.8) +
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = 0, ymax = 1), fill = pal[2], alpha = 0.8) +
  geom_rect(aes(xmin = x_max, xmax = max(metaData$Age) + 5, ymin = 0, ymax = 1), fill = pal[3], alpha = 0.8) +
  geom_vline(xintercept = x_min, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = x_max, linetype = "dashed", color = "black", size = 1) +
  geom_text(aes(x = min(metaData$Age) - 5 + (x_min - (min(metaData$Age) - 5))/2, y = 0.5, label = "Early-life"))+
  geom_text(aes(x = x_min + (x_max - x_min)/2, y = 0.5, label = "Mid-life"))+
  geom_text(aes(x = x_max + (max(metaData$Age) + 5 - x_max)/2, y = 0.5, label = "Late-life"))+
  theme_void()

# Make main plot
main <- ggplot() +
  geom_histogram(data = metaData, 
                 aes(x = Age, fill = AgeGroup), binwidth = 1, color = "grey") +
  #facet_grid(.~AgeGroup, scales = "free", space = "free") +
  geom_vline(xintercept = x_min, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = x_max, linetype = "dashed", color = "black", size = 1) +
  xlab("Age (years)") +
  ylab("Count") +
  xlim(c(min(metaData$Age) -5,max(metaData$Age) + 5)) +
  #ggtitle("Age Distribution") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"))

p <- ggarrange(top, main, nrow = 2, ncol = 1, heights = c(1,10), align = "v")

ggsave(p, file = "AgeDistribution.png", width = 8, height = 6)

