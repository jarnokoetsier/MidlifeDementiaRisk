# File: PreProcessing_ADNI.R
# Author: Jarno Koetsier
# Date: November 27, 2022

###############################################################################

# 0. Package and data management

###############################################################################

# Install Bioconductor packages
BiocManager::install(c("minfi", 
                       "wateRmelon", 
                       "IlluminaHumanMethylationEPICmanifest",
                       "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                       "minfiData",
                       "FlowSorted.Blood.EPIC"))

# Install CRAN packages
install.packages("RPMM")
install.packages("ggrepel")

# Install GitHub packages
devtools::install_github("markgene/maxprobes")

# Load packages
library(minfi)
library(wateRmelon)
library(tidyverse)
library(ggrepel)


# Clear workspace and console
rm(list = ls())
cat("\014") 

# Directories
DataDir <- "ADNI/"
OutputDir <- "ADNI/"

# Load data
load(paste0(DataDir,"RGset_ADNI.RData"))
load(paste0(DataDir,"RGset_test.RData"))
load(paste0(DataDir,"MetaData_ADNI.RData"))


RGset_all <- combineArrays(RGset, RGset_test, outType = "IlluminaHumanMethylationEPIC")
save(RGset_all,file = "RGset_ADNI_all.RData")
###############################################################################

# 1. Pre-normalization QC and filtering

###############################################################################


#=============================================================================#
# 1.1. Bisulfite conversion
#=============================================================================#

# Bisulfite conversion
bsc <- data.frame(bsc = wateRmelon::bscon(RGset_all))

# Keep samples with bisulfite conversion higher than 80%
keep <- rownames(bsc)[bsc$bsc > 80]
RGset_all <- RGset_all[,keep]

# Make bisulfite conversion plot
bisulfiteConv <- ggplot(bsc, aes(x = bsc)) +
  geom_histogram(bins = 30, color = "white") +
  theme_classic() +
  xlab("Median bisulfite conversion percentage") +
  ylab("Count") +
  ggtitle("Bisulfite Conversion") +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

ggsave(bisulfiteConv, file = paste0(OutputDir,"bisulfiteConv_ADNI.png"), height = 6, width = 8)

#=============================================================================#
# 1.2. Sample Quality
#=============================================================================#

# Get median intensities
qc <- getQC(preprocessRaw(RGset_all))
plotQC <- data.frame(ID = qc@rownames,
                     mMed = qc$mMed,
                     uMed = qc$uMed)

# Remove low quality samples
keep <- plotQC$ID[which(plotQC$uMed > -1*plotQC$mMed + 21)]
RGset_all <- RGset_all[,keep]

# Make plot
p <- ggplot() +
  geom_abline(intercept = 21, slope = -1, color = "grey", linewidth = 1.5, linetype = "dashed") +
  geom_point(data = plotQC, aes(x = mMed, y  = uMed, color = mMed*uMed), alpha = 0.5, size = 2) +
  #geom_text_repel(data = plotQC[setdiff(rownames(plotQC),keep),],aes(x = mMed, y = uMed, label = ID)) +
  xlab(expression("Median methylated "*log[2]*" intensity")) +
  ylab(expression("Median unmethylated "*log[2]*" intensity")) +
  ggtitle("Chipwide Median Intensity") +
  xlim(c(8,14)) +
  ylim(c(8,14)) +
  scale_color_viridis_c() +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

# Save plot
ggsave(p, file = paste0(OutputDir,"sampleQuality_ADNI.png"), width = 8, height = 6)



#=============================================================================#
# 1.3. Check sex labels
#=============================================================================#

# Predicted sex
predictedSex <- getSex(mapToGenome(RGset_all), cutoff = -2)
predictedSex <- data.frame(
  SampleID = predictedSex@rownames,
  PredictedSex = predictedSex$predictedSex,
  xMed = predictedSex$xMed,
  yMed = predictedSex$yMed
)
predictedSex <- inner_join(predictedSex,MetaData_baseline, by = c("SampleID" = "Basename"))
predictedSex$Sex <- ifelse(predictedSex$Sex == "M", "Male", "Female")

# Make plot
matchingSex <- ggplot(predictedSex) +
  geom_point(aes(x = xMed, y = yMed, color = Sex), alpha = 0.7) +
  geom_abline(intercept = -2, slope = 1, linetype = "dashed", linewidth = 1.5, color = "grey") +
  xlab(expression("Median X-chromosomal "*log[2]*" intensity")) +
  ylab(expression("Median Y-chromosomal "*log[2]*" intensity")) +
  ggtitle("Correctness of Sex Labels") +
  labs(color = "Sex label") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16)) +
  scale_color_brewer(palette = "Set1")

# Save plot
ggsave(matchingSex, file = paste0(OutputDir,"matchingSex_ADNI.png"), width = 8, height = 6)


###############################################################################

# 2. Normalization

###############################################################################

#single sample Noob (minfi)
set.seed(123)
methSet_noob1 <- preprocessNoob(RGset_all, 
                                offset = 15, 
                                dyeCorr = TRUE, 
                                verbose = FALSE,
                                dyeMethod = "single")

# BMIQ (wateRmelon)
set.seed(123)
methSet_allNorm <- wateRmelon::BMIQ(methSet_noob1, nfit = 10000)

# Save R objects
save(methSet_allNorm, file = paste0(DataDir,"methSet_allNorm_ADNI.RData"))


###############################################################################

# 3. Probe filtering

###############################################################################

gc()

# Load normalized data
load(paste0(DataDir,"methSet_allNorm_ADNI.RData"))

# Select probes that are in EXTEND
X_ADNI <- methSet_allNorm[rownames(X_test),]
sum(is.na(X_ADNI))
save(X_ADNI, file = paste0(DataDir,"X_ADNI.RData"))

# Get detection p-value
lumi_dpval <- detectionP(RGset_all, type = "m+u")
save(lumi_dpval, file = "lumi_dpval_ADNI.RData")










