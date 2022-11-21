# Install packages
BiocManager::install(c("minfi", 
                       "wateRmelon", 
                       "IlluminaHumanMethylationEPICmanifest",
                       "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))
install.packages("RPMM")
# Load packages
library(minfi)
library(wateRmelon)
library(tidyverse)

# Load data
load("rgsetAll.rdata")
load("metadat_extend.Rdata")

###############################################################################

# 1. Pre-processing and quality control

###############################################################################

#=============================================================================#
# 1.1. Filter samples
#=============================================================================#

# Check whether basename and ID are both unique identifiers
length(unique(dat$Basename))
length(unique(dat$ID))

# Get beta values from RGset object
bv <- getBeta(RGset)
length(intersect(RGset@colData@rownames,dat$Basename))

# Remove samples with age < 40 or age > 75 (keep midlife only)
keepSamples <- dat$Basename[(dat$Age >= 40) & (dat$Age <= 75)]

# Remove samples with missing values only
missingSample <- colnames(bv)[colSums(is.na(bv)) == nrow(bv)]
keepSamples <- setdiff(keepSamples, missingSample)

# Filter RGset
RGset1 <- RGset[,keepSamples]


#=============================================================================#
# 1.2. Pre-normalization quality control
#=============================================================================#

# Make a QC report
qcReport(RGset1, pdf= "qcReport.pdf")


# Outlier detection (wateRmelon)
outliers <- outlyx(RGset1, plot=TRUE)

# Bisulfite conversion
bsc <- data.frame(bsc = bscon(RGset1))
bisulfiteConv <- ggplot(bsc, aes(x = bsc)) +
  geom_histogram(bins = 30, color = "white") +
  theme_classic() +
  xlab("Median bisulfite conversion percentage") +
  ylab("Count") +
  ggtitle("Bisulfite Conversion") +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

ggsave(bisulfiteConv, file = "bisulfiteConv.png", height = 6, width = 8)

# Predicted sex
predictedSex <- getSex(mapToGenome(RGset1), cutoff = -2)
predictedSex <- data.frame(
  SampleID = predictedSex@rownames,
  PredictedSex = predictedSex$predictedSex,
  xMed = predictedSex$xMed,
  yMed = predictedSex$yMed
)
predictedSex <- inner_join(predictedSex,dat, by = c("SampleID" = "Basename"))
predictedSex$Sex <- ifelse(predictedSex$Sex == 1, "Male", "Female")

# Make plot
matchingSex <- ggplot(predictedSex) +
  geom_point(aes(x = xMed, y = yMed, color = Sex), alpha = 0.7) +
  geom_abline(intercept = -2, slope = 1, linetype = "dashed", linewidth = 1.5, color = "grey") +
  xlab(expression("Median X-chromosomal "*log[2]*" intensity")) +
  ylab(expression("Median Y-chromosomal "*log[2]*" intensity")) +
  labs(color = "Sex label") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

ggsave(matchingSex, file = "matchingSex.png", width = 8, height = 6)

#=============================================================================#
# 1.3. Normalization
#=============================================================================#

# 1) Remove probes with missing values
bv <- getBeta(RGset1)
RGset1 <- RGset1[rowSums(is.na(bv)) == 0,]


#single sample Noob (minfi) 
methSet_noob <- preprocessNoob(RGset1, 
                               offset = 15, 
                               dyeCorr = TRUE, 
                               verbose = FALSE,
                               dyeMethod = "single")

# BMIQ (wateRmelon)
methSet_allNorm <- BMIQ(methSet_noob)

#=============================================================================#
# 1.4. Remove probes
#=============================================================================#

# 1) Remove probes with high detection p-value:

# Get detection p-value
lumi_dpval <- detectionP(RGset1, type = "m+u")

# remove probes above detection p-value
lumi_failed <- lumi_dpval > 0.01

# Remove probes with detection p-value > 0.01 in at least one sample
# (Information leakage???)
lumi_dpval_keep <- rownames(lumi_failed)[rowSums(lumi_failed)==0]
methSet_allNorm1 <- methSet_allNorm[rownames(methSet_allNorm) %in% lumi_dpval_keep,]

# 2) Remove probes with SNPs
# drop the probes containing a SNP at the CpG interrogation and/or at the single nucleotide
# extension, for any minor allele frequency
snps <- getSnpInfo(RGset1)
keep <- snps@rownames[!(isTRUE(snps@listData$CpG_maf > 0) | isTRUE(snps@listData$SBE_maf > 0))]
methSet_allNorm1 <- methSet_allNorm1[rownames(methSet_allNorm1) %in% keep,]


# 4) Remove probes that map to X and Y chromosome
annotation <- getAnnotation(RGset1)
keep <- annotation@rownames[(annotation@listData$chr != "chrX") & (annotation@listData$chr != "chrY")]
methSet_allNorm1 <- methSet_allNorm1[rownames(methSet_allNorm1) %in% keep,]


# Check success of normalization
plotDensity <- gather(as.data.frame(methSet_allNorm1))
densityNorm <- ggplot(plotDensity) +
  geom_density(aes(x = value, color = key)) +
  xlab("Beta value") +
  ylab("Density") +
  ggtitle("Post-normalization") +
  scale_color_viridis_c() +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16)) 
ggsave(densityNorm, file = "densityNorm.png", width = 8, height = 6)

plotDensity <- gather(as.data.frame(bv))
densityRaw <- ggplot(plotDensity) +
  geom_density(aes(x = value, color = key)) +
  xlab("Beta value") +
  ylab("Density") +
  ggtitle("Pre-normalization") +
  scale_color_viridis_c() +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))
ggsave(densityRaw, file = "densityRaw.png", width = 8, height = 6)


