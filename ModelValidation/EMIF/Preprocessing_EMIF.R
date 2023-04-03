# Load packages
library(minfi)
library(wateRmelon)
library(tidyverse)
library(ggrepel)
library(RPMM)


# Clear workspace and console
rm(list = ls())
cat("\014") 

# Directories
DataDir <- "EMIF/"
OutputDir <- "EMIF/"

# Load data
load(paste0(DataDir,"RGset_EMIF1.RData"))
load(paste0(DataDir,"RGset_EMIF2.RData"))
load(paste0(DataDir,"RGset_EMIF3.RData"))
load(paste0(DataDir,"RGset_EMIF4.RData"))

RGset_all <- combineArrays(RGset_EMIF1, RGset_EMIF2,
                           outType = "IlluminaHumanMethylationEPIC")
RGset_all1 <- combineArrays(RGset_EMIF3, RGset_EMIF4,
                           outType = "IlluminaHumanMethylationEPIC")

save(RGset_all, RGset_all1,file = "EMIF/RGset_EMIFF_all.RData")

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

bsc <- data.frame(bsc = wateRmelon::bscon(RGset_all1))

# Keep samples with bisulfite conversion higher than 80%
keep <- rownames(bsc)[bsc$bsc > 80]
RGset_all1 <- RGset_all1[,keep]




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

qc <- getQC(preprocessRaw(RGset_all1))
plotQC <- data.frame(ID = qc@rownames,
                     mMed = qc$mMed,
                     uMed = qc$uMed)

# Remove low quality samples
keep <- plotQC$ID[which(plotQC$uMed > -1*plotQC$mMed + 21)]
RGset_all1 <- RGset_all1[,keep]

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
methSet_noob <- preprocessNoob(RGset_all, 
                                offset = 15, 
                                dyeCorr = TRUE, 
                                verbose = FALSE,
                                dyeMethod = "single")

# BMIQ (wateRmelon)
set.seed(123)
methSet_allNorm <- wateRmelon::BMIQ(methSet_noob, nfit = 10000)

# Save R objects
save(methSet_allNorm, file = paste0(DataDir,"methSet_allNorm_EMIF.RData"))

#single sample Noob (minfi)
set.seed(123)
methSet_noob1 <- preprocessNoob(RGset_all1, 
                                offset = 15, 
                                dyeCorr = TRUE, 
                                verbose = FALSE,
                                dyeMethod = "single")

# BMIQ (wateRmelon)
set.seed(123)
methSet_allNorm1 <- wateRmelon::BMIQ(methSet_noob1, nfit = 10000)

# Save R objects
save(methSet_allNorm1, file = paste0(DataDir,"methSet_allNorm_EMIF1.RData"))


###############################################################################

# 3. Probe filtering

###############################################################################

gc()

# Load normalized data
load(paste0(DataDir,"methSet_allNorm_EMIF.RData"))
load(paste0(DataDir,"methSet_allNorm_EMIF1.RData"))

X_all <- cbind(methSet_allNorm, methSet_allNorm1)

# Select probes that are in EXTEND
load("~/Data/X_test.RData")
X_EMIF <- X_all[rownames(X_test),]
sum(is.na(X_EMIF))
save(X_EMIF, file = paste0(DataDir,"X_EMIF.RData"))


# Get detection p-value
lumi_dpval1 <- detectionP(RGset_all, type = "m+u")
lumi_dpval2 <- detectionP(RGset_all1, type = "m+u")
lumi_dpval <- cbind(lumi_dpval1, lumi_dpval2)
save(lumi_dpval, file = paste0(DataDir,"lumi_dpval_EMIF.RData"))

#PCA
pcaList <-  prcomp(t(X_EMIF),        
                   retx = TRUE,
                   center =TRUE,
                   scale = TRUE,
                   rank. = 10)

# Save pcaList object
save(pcaList, file = "EMIF/pcaList_EMIF.RData")

# Get PCA scores
PCAscores <- as.data.frame(pcaList$x)
PCAscores$ID <- rownames(PCAscores)

# Get explained variance
explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)

# Make PCA score: PC1 vs PC2
p_12 <- ggplot(data = PCAscores, aes(x = PC1, y = PC2)) +
  geom_point() +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text())

ggsave(p_12, file = "EMIF/PCA_EMIF.png", width = 8, height = 6)

keepSamples <- rownames(PCAscores)[PCAscores$PC1 > -2000]
X_EMIF <- X_EMIF[,keepSamples]
lumi_dpval <- lumi_dpval[,keepSamples]
  
metaData_EMIF <- read_csv("EMIF/EMIF_metadata.csv")
save(metaData_EMIF, file = "EMIF/metaData_EMIF.RData")
samples <- intersect(colnames(X_EMIF), metaData_EMIF$X)

