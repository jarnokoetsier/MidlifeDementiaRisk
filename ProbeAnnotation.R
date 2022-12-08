
# Load packages
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(tidyverse)

# Get probe annotations
ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annList <- ann@listData
probe_annotation <- data.frame(
  ID = annList$Name,
  Chr = annList$chr,
  Position = annList$pos,
  Strand = annList$strand,
  Type = annList$Type,
  CpG_MAF = annList$CpG_maf,
  Relation_to_Island = annList$Relation_to_Island,
  Gene_Name = annList$UCSC_RefGene_Name,
  Gene_ID = annList$UCSC_RefGene_Accession,
  Gene_Group = annList$UCSC_RefGene_Group,
  Regulatory_Group = annList$Regulatory_Feature_Group
)

# Divide probes into three classes:
# 1. Promotor
# 2. Gene body
# 3. Intergenic

body <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "Body")]
tss200 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "TSS200")]
tss1500 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "TSS1500")]
utr5 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "5'UTR")]
utr3 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "3'UTR")]
exon1 <- probe_annotation$ID[str_detect(probe_annotation$Gene_Group, "1stExon")]

promotor <- unique(c(tss200,tss1500))
geneBody <- setdiff(unique(c(utr5,exon1, body, utr3)),promotor)
interGenic <- setdiff(probe_annotation$ID, unique(c(promotor, geneBody)))

probe_annotation$Class <- NA
probe_annotation$Class[probe_annotation$ID %in% promotor] <- "Promotor"
probe_annotation$Class[probe_annotation$ID %in% geneBody] <- "Gene Body"
probe_annotation$Class[probe_annotation$ID %in% interGenic] <- "Intergenic"

# Save as .RData object
save(probe_annotation,file = "E:/Thesis/MLData/probe_annotation.RData")
