# Install packages

# Bioconductor packages
BiocManager::install(c("minfi", 
                       "wateRmelon", 
                       "IlluminaHumanMethylationEPICmanifest",
                       "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                       "minfiData",
                       "FlowSorted.Blood.EPIC"))

# CRAN packages
install.packages("RPMM")

# GitHub packages
devtools::install_github("markgene/maxprobes")

# Load packages
library(minfi)
library(wateRmelon)
library(tidyverse)
library(maxprobes)

# Load data
load("rgsetAll.rdata")
load("metadat_extend.Rdata")

###############################################################################

# 1. Pre-normalization QC and filtering

###############################################################################

# Check whether basename and ID are both unique identifiers
length(unique(dat$Basename))
length(unique(dat$ID))

# Get beta values from RGset object
bv <- getBeta(RGset)
length(intersect(RGset@colData@rownames,dat$Basename))

#=============================================================================#
# 1.1. Age filtering
#=============================================================================#

# Remove samples with age < 40 or age > 75 (keep midlife only)
keepSamples <- dat$Basename[(dat$Age >= 40) & (dat$Age <= 75)]

# Remove samples with missing values only
missingSample <- colnames(bv)[colSums(is.na(bv)) == nrow(bv)]
keepSamples <- setdiff(keepSamples, missingSample)

# Filter RGset
RGset1 <- RGset[,keepSamples]


# Make a QC report (NOT NECESARY TO RUN)
qcReport(RGset1, pdf= "qcReport.pdf")


# Outlier detection (wateRmelon) (NOT NECESARY TO RUN)
outliers <- outlyx(RGset1, plot=TRUE)


#=============================================================================#
# 1.2. Bisulfite conversion
#=============================================================================#

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

#=============================================================================#
# 1.3. Sample Quality
#=============================================================================#

# Get median intensities
qc <- getQC(preprocessRaw(RGset1))
plotQC <- data.frame(ID = qc@rownames,
                     mMed = qc$mMed,
                     uMed = qc$uMed)

# Make plot
p <- ggplot(plotQC) +
  geom_point(aes(x = mMed, y  = uMed, color = mMed*uMed), alpha = 0.5, size = 2) +
  geom_abline(intercept = 21, slope = -1, color = "grey", linewidth = 1.5, linetype = "dashed") +
  xlab(expression("Median methylated "*log[2]*" intensity")) +
  ylab(expression("Median unmethylated "*log[2]*" intensity")) +
  xlim(c(8,14)) +
  ylim(c(8,14)) +
  scale_color_viridis_c() +
  theme_minimal() +
  theme(legend.position = "none")

# Save plot
ggsave(p, file = "sampleQuality.png", width = 8, height = 6)


#=============================================================================#
# 1.4. Check sex labels
#=============================================================================#

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

# Save plot
ggsave(matchingSex, file = "matchingSex.png", width = 8, height = 6)


#=============================================================================#
# 1.5. Remove probes with missing values
#=============================================================================#

# Remove probes with missing values
bv <- getBeta(RGset1)
RGset1 <- RGset1[rowSums(is.na(bv)) == 0,]



###############################################################################

# 2. Normalization

###############################################################################


#single sample Noob (minfi) 
methSet_noob <- preprocessNoob(RGset1, 
                               offset = 15, 
                               dyeCorr = TRUE, 
                               verbose = FALSE,
                               dyeMethod = "single")

# BMIQ (wateRmelon) (nfit = 50000)
methSet_allNorm <- BMIQ(methSet_noob, nfit = 50000)

# Save R objects
save(methSet_allNorm, file = "methSet_allNorm.RData")
save(RGset1, file = "RGset1.RData")

###############################################################################

# 3. Post-normalization QC and filtering

###############################################################################

#=============================================================================#
# 3.1. Detection P-value
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


#=============================================================================#
# 3.2. SNPs
#=============================================================================#

# 2) Remove probes with SNPs
# drop the probes containing a SNP at the CpG interrogation and/or at the single nucleotide
# extension, for any minor allele frequency
snps <- getSnpInfo(RGset1)
keep <- snps@rownames[!(isTRUE(snps@listData$CpG_maf > 0) | isTRUE(snps@listData$SBE_maf > 0))]
methSet_allNorm1 <- methSet_allNorm1[rownames(methSet_allNorm1) %in% keep,]


#=============================================================================#
# 3.3.  Sex chromosomes
#=============================================================================#

# 3) Remove probes that map to X and Y chromosome
annotation <- getAnnotation(RGset1)
keep <- annotation@rownames[(annotation@listData$chr != "chrX") & (annotation@listData$chr != "chrY")]
methSet_allNorm1 <- methSet_allNorm1[rownames(methSet_allNorm1) %in% keep,]


#=============================================================================#
# 3.4. Cross-reactive probes
#=============================================================================#

# Remove cross reactive probes
remove <- maxprobes::xreactive_probes(array_type = "EPIC")
keep <- setdiff(rownames(methSet_allNorm1), remove)
methSet_allNorm1 <- methSet_allNorm1[rownames(methSet_allNorm1) %in% keep,]
save(methSet_allNorm1, file = "methSet_allNorm1.RData")


#=============================================================================#
# 3.5. Check normalization 
#=============================================================================#

# Plot density after normalization
plotDensity <- gather(as.data.frame(methSet_allNorm1))
densityNorm <- ggplot(plotDensity) +
  geom_density(aes(x = value, color = key)) +
  xlab("Beta value") +
  ylab("Density") +
  ggtitle("Post-normalization") +
  scale_color_viridis_d() +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16)) 
# Save plot
ggsave(densityNorm, file = "densityNorm.png", width = 8, height = 6)

# Plot density before normalization
plotDensity <- gather(as.data.frame(bv))
densityRaw <- ggplot(plotDensity) +
  geom_density(aes(x = value, color = key)) +
  xlab("Beta value") +
  ylab("Density") +
  ggtitle("Pre-normalization") +
  scale_color_viridis_d() +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))
# Save plot
ggsave(densityRaw, file = "densityRaw.png", width = 8, height = 6)


#=============================================================================#
# 3.5. Check for outliers
#=============================================================================#

# 1. PCA

# Perform PCA
pcaList <-  prcomp(t(methSet_allNorm1),        
                   retx = TRUE,
                   center =TRUE,
                   scale = TRUE,
                   rank. = 10)

save(pcaList, file = "pcaList.RData")

# Get PCA scores
PCAscores <- as.data.frame(pcaList$x)
PCAscores$ID <- rownames(PCAscores)
PCAscores <- inner_join(PCAscores, dat, by = c("ID" = "Basename"))

# Get explained variance
explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)

# Make PCA score
p <- ggplot(data = PCAscores, aes(x = PC1, y = PC2, color = Age)) +
  geom_point(alpha = 0.9, size = 2) +
  ggtitle(label = "Age") +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_color_viridis_c()

# save plot
ggsave(p, file = "PCAplot.png", width = 8, height = 6)


# 2. Distance-Distance plot

# Split data into reconstruction and residuals based on the first ten PC's
loadings <- pcaList$rotation
reconstruction <- as.matrix(PCAscores[,1:10]) %*% t(loadings[,1:10])
residuals <- t(methSet_allNorm1) - reconstruction

# Calculate the orthogonal distances
ortDist<- sqrt(rowSums(residuals^2))

# Calculate the mahalanobis distances
coVar <- cov(PCAscores[,1:10])
center <- colMeans(PCAscores[,1:10])
mahalanobisDist <- mahalanobis(as.matrix(PCAscores[,1:10]), center = center, cov = coVar)

# Save distance metrics into a data frame
distanceDF <- data.frame(OD = ortDist,
                         MD = mahalanobisDist,
                         Age = PCAscores$Age)

# Distance-distance plot, colored by age
DDplot <- ggplot(distanceDF, aes(x = MD, y = OD, color = Age)) +
  geom_point(alpha = 0.9, size = 2) +
  xlab("Score distance (10 PC)") + 
  ylab("Orthogonal distance") +
  labs(color = "Age") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "right") +
  scale_color_viridis_c()

# Save plot
ggsave(DDplot, file = "DDplot.png", width = 8, height = 6)


# 3. Explained variances

# Put explained variances into a data frame
explVarDF <- data.frame(PC = factor(paste0("PC",1:10), levels = paste0("PC",1:10)),
                        Variance = explVar[1:10])

# Make plot
p <- ggplot(explVarDF) +
  geom_bar(aes(x = PC, y = Variance, fill = Variance), 
           alpha = 0.8, color = "black", stat = "identity") +
  xlab(NULL) +
  ylab("Explained Variance (%)") +
  scale_fill_viridis_c() +
  theme_classic() +
  theme(legend.position = "none")

# Save plot
ggsave(p, file = "explVar.png", width = 8, height = 6)
  
###############################################################################

# 4. Estimate cell type composition

###############################################################################

# Estimate cell type composition
library(FlowSorted.Blood.EPIC)

cellType <-estimateCellCounts.wmln(as.methylumi(preprocessRaw(RGset1)),
                                   platform = "EPIC",
                                   mn = NULL,
                                   un = NULL,
                                   bn = NULL,
                                   perc = 1,
                                   compositeCellType = "Blood",
                                   probeSelect = "auto",
                                   cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
                                   referencePlatform = "IlluminaHumanMethylationEPIC",
                                   returnAll = FALSE,
                                   meanPlot = FALSE,
                                   verbose= FALSE)

cellType <- as.data.frame(cellType)
cellType$ID <- rownames(cellType)
save(cellType, file = "cellType.RData")
PCAscores <- inner_join(PCAscores, cellType, by = c("ID" = "ID"))

# Get explained variance
explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)

# Make PCA score
p <- ggplot(data = PCAscores, aes(x = PC1, y = PC2, color = CD8T)) +
  geom_point(alpha = 0.9, size = 2) +
  ggtitle(label = "CD8T") +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_blank()) +
  scale_color_viridis_c()

# save plot
ggsave(p, file = "PCAplot_CD8T.png", width = 8, height = 6)
