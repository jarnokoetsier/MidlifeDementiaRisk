# Load packages
library(tidyverse)
library(data.table)
library(R.utils)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37) # Bioconductor package
library(colochelpR) # GitHub package
library(rtracklayer) # Bioconductor package


# Clear workspace and console
rm(list = ls())
cat("\014") 

#=============================================================================#
# FILL IN:

# set directories
rawDir <- "E:/Thesis/GWAS/Raw/"
curDir <- "E:/Thesis/GWAS/Curated/"

# File name of raw file
fileName <- "SBPmanual_UKB.tsv"

# Information
trait <- "SBPmanual"
year <- 2018
author <- "UKBB"
#=============================================================================#

# Read raw file
dataObj <- fread(paste0(rawDir,fileName))

# Column names:
# Continuous:
allNames_con <- c("CHR", "BP", "A1", "A2", "N", "BETA", "P")

# Categorical:
allNames_cat <- c("CHR", "BP", "A1", "A2", "N", "OR", "P")

# Capital letters for SNPs (if needed)
dataObj$A1 <- toupper(dataObj$A1)
dataObj$A2 <- toupper(dataObj$A2)

# Add number of samples (if needed)
dataObj$N <- rep(807553,nrow(dataObj))

# Convert rs id to Chr:Bp
dataObj <- convert_rs_to_loc(dataObj, "MarkerName", dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh37)

# Split chromosome and base pair (if needed)
dataObj$CHR <- str_remove_all(unlist(lapply(str_split(dataObj$variant, ":"), `[[`, 1)), "chr")
dataObj$BP <- unlist(lapply(str_split(dataObj$variant, ":"), `[[`, 2))
dataObj$A1 <- unlist(lapply(str_split(dataObj$variant, ":"), `[[`, 3))
dataObj$A2 <- unlist(lapply(str_split(dataObj$variant, ":"), `[[`, 4))

# Extract columns
dataObj_cur <- dataObj[, c("CHR", "BP", "A1", "A2", "N", "beta", "pval")]
dataObj_cur <- dataObj[, c("CHR", "BP", "A1", "A2", "N", "OR", "P")]

# Convert beta to OR (only for categorical data)
dataObj_cur$BETA <- exp(dataObj_cur$BETA)

# correct col names
colnames(dataObj_cur) <- allNames_con

#=============================================================================#
# Convert genome build (if needed)

# Find chain file at the liftOver tab at the UCSC website:
# http://hgdownload.soe.ucsc.edu/downloads.html

# import the chain file
chainObject <- import.chain("GWAS/hg18ToHg19.over.chain")

# specify coordinates to liftover
dataObj <- dataObj[!is.na(dataObj$BP),]
grObject <- GRanges(seqnames=paste0("chr", dataObj$CHR), 
                    ranges=IRanges(start=dataObj$BP, end=dataObj$BP))


# run liftOver
results <- as.data.frame(liftOver(grObject, chainObject))

# Convert genomic location
dataObj_cur <- dataObj[results$group,]
dataObj_cur$CHR<- str_remove_all(results$seqnames, "chr")
dataObj_cur$BP <- results$start
#=============================================================================#

# Write output
write.table(dataObj_cur, file = paste0(curDir, trait, "_", year, "_", author, "_curated.summaries"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

