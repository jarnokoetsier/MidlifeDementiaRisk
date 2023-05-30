###############################################################################

# Estimate cell type composition

###############################################################################
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
install.packages("ggrepel")

# GitHub packages
devtools::install_github("markgene/maxprobes")

# Load packages
library(minfi)
library(wateRmelon)
library(tidyverse)
library(maxprobes)
library(ggrepel)
library(FlowSorted.Blood.EPIC)

# Load data
load("metaData_ageFil.RData")
load("rgsetAll.rdata")
RGset <- RGset[,dat$Basename]

# Estimate cell type composition
cellType1 <- estimateCellCounts.wmln(as.methylumi(preprocessRaw(RGset)),
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

# Format data
cellType <- as.data.frame(cellType1)
cellType$ID <- rownames(cellType)

# Save data
save(cellType, file = "cellType.RData")

