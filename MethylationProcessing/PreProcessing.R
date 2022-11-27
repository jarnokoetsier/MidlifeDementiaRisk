# File: PreProcessing.R
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
library(maxprobes)
library(ggrepel)

# Load data
load("rgsetAll.rdata")
load("metaData_ageFil.RData")

###############################################################################

# 1. Pre-normalization QC and filtering

###############################################################################

# Check whether basename and ID are both unique identifiers
length(unique(dat$Basename))
length(unique(dat$ID))

# Get beta values from RGset object
RGset <- RGset[,dat$Basename]
bv <- getBeta(RGset)

#=============================================================================#
# 1.1. Sample filtering
#=============================================================================#

# Remove samples with more than 20% missing missing values 
keepSamples <- colnames(bv)[colSums(is.na(bv)) < 0.2*nrow(bv)]

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
  ggtitle("Correctness of Sex Labels") +
  labs(color = "Sex label") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16)) +
  scale_color_brewer(palette = "Set1")

# Save plot
ggsave(matchingSex, file = "matchingSex.png", width = 8, height = 6)


###############################################################################

# 2. Normalization

###############################################################################
load("RGset1.RData")

#single sample Noob (minfi) 
methSet_noob1 <- preprocessNoob(RGset1, 
                               offset = 15, 
                               dyeCorr = TRUE, 
                               verbose = FALSE,
                               dyeMethod = "single")

# BMIQ (wateRmelon) (nfit = 50000)
methSet_allNorm <- BMIQ(methSet_noob1, nfit = 10000)

# Save R objects
save(methSet_allNorm, file = "methSet_allNorm1.RData")

###############################################################################

# 3. Post-normalization probe filtering

###############################################################################

load("methSet_allNorm1.RData")

#=============================================================================#
# 3.1. Detection P-value
#=============================================================================#

# 1) Remove probes with high detection p-value:

# Get detection p-value
lumi_dpval <- detectionP(RGset1, type = "m+u")

# remove probes above detection p-value
lumi_failed <- lumi_dpval > 0.01

# Remove probes with detection p-value > 0.01 in at least one sample
lumi_dpval_keep <- rownames(lumi_failed)[rowSums(lumi_failed)==0]
methSet_allNorm_fil <- methSet_allNorm[rownames(methSet_allNorm) %in% lumi_dpval_keep,]


#=============================================================================#
# 3.2. SNPs
#=============================================================================#

# 2) Remove probes with SNPs
# drop the probes containing a SNP at the CpG interrogation and/or at the single nucleotide
# extension, for any minor allele frequency
snps <- getSnpInfo(RGset1)
keep <- snps@rownames[!(isTRUE(snps@listData$CpG_maf > 0) | isTRUE(snps@listData$SBE_maf > 0))]
methSet_allNorm_fil <- methSet_allNorm_fil[rownames(methSet_allNorm_fil) %in% keep,]


#=============================================================================#
# 3.3. Cross-reactive probes
#=============================================================================#

# Remove cross reactive probes
remove <- maxprobes::xreactive_probes(array_type = "EPIC")
keep <- setdiff(rownames(methSet_allNorm_fil), remove)
methSet_allNorm_fil <- methSet_allNorm_fil[rownames(methSet_allNorm_fil) %in% keep,]


#=============================================================================#
# 3.4. Missing/unrealistic values & non-unique rows
#=============================================================================#

# Check for NA values
sum(is.na(methSet_allNorm_fil))

# Check if values are lower or equal to zero
sum(methSet_allNorm_fil <= 0)

# Check if values are larger or equal to one
sum(methSet_allNorm_fil >= 1)

# Check for non-unique probes
length(unique(rownames(methSet_allNorm_fil))) == nrow(methSet_allNorm_fil)
nrow(methSet_allNorm_fil)

# Save methylation matrix
save(methSet_allNorm_fil, file = "methSet_allNorm_fil.RData")


###############################################################################

# 4. Quality control

###############################################################################

load("methSet_allNorm_fil.RData")
load("cellType.RData")

#=============================================================================#
# 4.1. Density plots
#=============================================================================#

# Plot density after normalization
plotDensity <- gather(as.data.frame(methSet_allNorm_fil))
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
# 4.2. PCA plot
#=============================================================================#

# Perform PCA
pcaList <-  prcomp(t(methSet_allNorm_fil),        
                   retx = TRUE,
                   center =TRUE,
                   scale = TRUE,
                   rank. = 10)

# Save pcaList object
save(pcaList, file = "pcaList.RData")

# Get PCA scores
PCAscores <- as.data.frame(pcaList$x)
PCAscores$ID <- rownames(PCAscores)
PCAscores <- inner_join(PCAscores, dat, by = c("ID" = "Basename"))

# Combine with cell type composition
PCAscores <- inner_join(PCAscores,cellType, by = c("ID" = "ID"))

# Get explained variance
explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)

# Make PCA score: PC1 vs PC2
p_12 <- ggplot(data = PCAscores, aes(x = PC1, y = PC2)) +
  stat_ellipse(geom = "polygon",
               fill = "red",
               type = "norm", 
               alpha = 0.25,
               level = 0.95) +
  stat_ellipse(geom = "polygon",
               color = "red",
               alpha = 0,
               linetype = 2,
               type = "norm",
               level = 0.99) +
  geom_point(alpha = 0.9, size = 2, aes(color = CD8T)) +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  labs(color = "CD8 T-cells") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text()) +
  scale_color_viridis_c()

# save plot
ggsave(p_12, file = "PCAplot_1vs2.png", width = 8, height = 6)

# Make PCA score plot: PC3 vs PC4
p_34 <- ggplot(data = PCAscores, aes(x = PC3, y = PC4)) +
  stat_ellipse(geom = "polygon",
               fill = "red",
               type = "norm", 
               alpha = 0.25,
               level = 0.95) +
  stat_ellipse(geom = "polygon",
               color = "red",
               alpha = 0,
               linetype = 2,
               type = "norm",
               level = 0.99) +
  geom_point(alpha = 0.9, size = 2, aes(color = CD4T)) +
  xlab(paste0("PC3 (", explVar[3],"%)")) +
  ylab(paste0("PC4 (", explVar[4],"%)")) +
  labs(color = "CD4 T-cells") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text()) +
  scale_color_viridis_c()

# save plot
ggsave(p_34, file = "PCAplot_3vs4.png", width = 8, height = 6)

# Make PCA score plot: PC5 vs PC6
PCAscores$Sex <- ifelse(PCAscores$Sex == 1, "Male", "Female")
p_56 <- ggplot(data = PCAscores, aes(x = PC5, y = PC6)) +
  stat_ellipse(geom = "polygon",
               fill = "red",
               type = "norm", 
               alpha = 0.25,
               level = 0.95) +
  stat_ellipse(geom = "polygon",
               color = "red",
               alpha = 0,
               linetype = 2,
               type = "norm",
               level = 0.99) +
  geom_point(alpha = 0.9, size = 2, aes(color = Sex)) +
  xlab(paste0("PC5 (", explVar[5],"%)")) +
  ylab(paste0("PC6 (", explVar[6],"%)")) +
  labs(color = "Sex") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text()) +
  scale_color_brewer(palette = 'Set1')

# save plot
ggsave(p_56, file = "PCAplot_5vs6.png", width = 8, height = 6)


#=============================================================================#
# 4.3. Distance-distance plot
#=============================================================================#

# Split data into reconstruction and residuals based on the first ten PC's
loadings <- pcaList$rotation
reconstruction <- as.matrix(PCAscores[,1:10]) %*% t(loadings[,1:10])
residuals <- t(methSet_allNorm_fil) - reconstruction

# Calculate the orthogonal distances
ortDist<- sqrt(rowSums(residuals^2))

# Calculate the mahalanobis distances
coVar <- cov(t(t(PCAscores[,1:10])*(explVar[1:10]/100)))
center <- colMeans(t(t(PCAscores[,1:10])*(explVar[1:10]/100)))
mahalanobisDist <- mahalanobis(as.matrix(PCAscores[,1:10]), center = center, cov = coVar)

# Save distance metrics into a data frame
distanceDF <- data.frame(OD = ortDist,
                         MD = mahalanobisDist,
                         Age = PCAscores$Age,
                         ID = PCAscores$ID)

# Distance-distance plot, colored by age
DDplot <- ggplot() +
  geom_point(data = distanceDF, aes(x = MD, y = OD, color = Age), 
             alpha = 0.9, size = 2) +
  geom_text_repel(data = distanceDF[(distanceDF$MD > 690000) |
                                      (distanceDF$OD > 2000),], 
            aes(x = MD, y = OD, label = ID), size = 2) +
  #geom_vline(xintercept = log(500000), linetype = "dashed", color = 'red') +
  #geom_hline(yintercept = log(1500), linetype = "dashed", color = "red") +
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

# Label the identified "outliers" in regular PCA score plot
outliers <- distanceDF$ID[(distanceDF$MD > 690000) | (distanceDF$OD > 2000)]

# PC1 vs PC2
p_12_label <- ggplot() +
  stat_ellipse(geom = "polygon",
               data = PCAscores,
               aes(x = PC1, y = PC2),
               fill = "red",
               type = "norm", 
               alpha = 0.25,
               level = 0.95) +
  stat_ellipse(geom = "polygon",
               data = PCAscores,
               aes(x = PC1, y = PC2),
               color = "red",
               alpha = 0,
               linetype = 2,
               type = "norm",
               level = 0.99) +
  geom_point(data = PCAscores, aes(x = PC1, y = PC2, color = CD8T), 
             alpha = 0.9, size = 2) +
  geom_text_repel(data = PCAscores[PCAscores$ID %in% outliers,], 
                  aes(x = PC1, y = PC2, label = ID), size = 2) +
  xlab(paste0("PC1 (", explVar[1],"%)")) +
  ylab(paste0("PC2 (", explVar[2],"%)")) +
  labs(color = "CD8 T-cells") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text()) +
  scale_color_viridis_c()

# save plot
ggsave(p_12_label, file = "PCAplot_1vs2_label.png", width = 8, height = 6)

# PC3 vs PC4
p_34_label <- ggplot() +
  stat_ellipse(geom = "polygon",
               data = PCAscores,
               aes(x = PC3, y = PC4),
               fill = "red",
               type = "norm", 
               alpha = 0.25,
               level = 0.95) +
  stat_ellipse(geom = "polygon",
               data = PCAscores,
               aes(x = PC3, y = PC4),
               color = "red",
               alpha = 0,
               linetype = 2,
               type = "norm",
               level = 0.99) +
  geom_point(data = PCAscores, aes(x = PC3, y = PC4, color = CD4T), 
             alpha = 0.9, size = 2) +
  geom_text_repel(data = PCAscores[PCAscores$ID %in% outliers,], 
                  aes(x = PC3, y = PC4, label = ID), size = 2) +
  xlab(paste0("PC3 (", explVar[3],"%)")) +
  ylab(paste0("PC4 (", explVar[4],"%)")) +
  labs(color = "CD4 T-cells") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text()) +
  scale_color_viridis_c()

# save plot
ggsave(p_34_label, file = "PCAplot_3vs4_label.png", width = 8, height = 6)

# PC5 vs PC6
p_56_label <- ggplot() +
  stat_ellipse(geom = "polygon",
               data = PCAscores,
               aes(x = PC5, y = PC6),
               fill = "red",
               type = "norm", 
               alpha = 0.25,
               level = 0.95) +
  stat_ellipse(geom = "polygon",
               data = PCAscores,
               aes(x = PC5, y = PC6),
               color = "red",
               alpha = 0,
               linetype = 2,
               type = "norm",
               level = 0.99) +
  geom_point(data = PCAscores, aes(x = PC5, y = PC6, color = Sex), 
             alpha = 0.9, size = 2) +
  geom_text_repel(data = PCAscores[PCAscores$ID %in% outliers,], 
                  aes(x = PC5, y = PC6, label = ID), size = 2) +
  xlab(paste0("PC5 (", explVar[5],"%)")) +
  ylab(paste0("PC6 (", explVar[6],"%)")) +
  labs(color = "Sex") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10),
        legend.position = "bottom",
        legend.title = element_text()) +
  scale_color_brewer(palette = "Set1")

# save plot
ggsave(p_56_label, file = "PCAplot_5vs6_label.png", width = 8, height = 6)


#=============================================================================#
# 4.4. PCA explained variances
#=============================================================================#

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
  

#=============================================================================#
# 4.5. PCA correlations
#=============================================================================#

PCs <- PCAscores[,1:10]
meta <- PCAscores[,c("Age", "Sex", colnames(cellType))]
meta <- meta[,-c(9,10)]
meta$Sex <- ifelse(meta$Sex == "Male", 0,1)
colnames(meta) <- c("Age", "Sex", "CD8 T-cells", "CD4 T-cells", "NK-cells", "B-cells", "Monocytes",
                    "Neutrophils")

corDF <- expand.grid(colnames(PCs), colnames(meta))
colnames(corDF) <- c("PCs", "Meta")
corDF$Correlation <- rep(NA, nrow(corDF))
for (i in 1:nrow(corDF)){
  corDF$Correlation[i] <- cor(PCs[,corDF$PCs[i]], meta[,corDF$Meta[i]], method = "spearman")
}

p <- ggplot() +
  geom_point(data = corDF, aes(x = PCs, y = Meta, color = Correlation, size = abs(Correlation))) +
  labs(color = "Spearman\nCorrelation", size = "|Spearman\nCorrelation|") +
  scale_color_viridis_c(limits = c(-1, 1), oob = scales::squish) +
  theme(axis.title = element_blank())

ggsave(p,file = "correlationPlot.png", width = 8, height = 6)

