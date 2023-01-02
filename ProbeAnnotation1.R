# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(tidyverse)
library(foreach)
library(doParallel)

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load probe annotation data
load("probe_annotation.RData")

# Filter for probes
probe_annotation <- probe_annotation[(probe_annotation$Chr != "chrY") & (probe_annotation$Chr != "ChrX"),]


# Get associated gene(s)

Probe2Gene <- NULL
for (i in 1:nrow(probe_annotation)){
  temp <- data.frame(
    Gene = unlist(str_split(probe_annotation$Gene_Name[i], ";")),
    CpG = rep(probe_annotation$ID[i],length(unlist(str_split(probe_annotation$Gene_Name[i], ";")))),
    Association = unlist(str_split(probe_annotation$Gene_Group[i], ";"))
  )
  Probe2Gene <- rbind.data.frame(Probe2Gene,temp)
}


# Register cores for parallel computing
nCores <- 3
cl <- makeCluster(nCores)
registerDoParallel(cl)
Probe2Gene <- foreach (i = 1:nrow(probe_annotation), 
                       .packages = "tidyverse", 
                       .combine = rbind.data.frame, 
                       .inorder = FALSE) %dopar% {
  temp <- data.frame(
    Gene = unlist(str_split(probe_annotation$Gene_Name[i], ";")),
    CpG = rep(probe_annotation$ID[i],length(unlist(str_split(probe_annotation$Gene_Name[i], ";")))),
    Association = unlist(str_split(probe_annotation$Gene_Group[i], ";"))
  )
  return(temp)
}
# Stop clusters
stopCluster(cl)