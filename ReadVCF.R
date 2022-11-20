install.packages("vcfR")

# Load vcfR package
library(vcfR)

# Read vcf file
vcf <- read.vcfR(file = "chr19.dose.vcf.gz")

# Investigate structure of vcf object
head(vcf@gt[,1:5])
head(vcf@fix)

# Get positions of SNPs
positions <- vcf@fix

# Get meta data 
meta <- vcf@meta

# hg19, GRCh37
#19:45411941:T:C (rs429358)
#19:45412079:C:T (rs7412)

# Positions relating to APOE4 status
positions[positions[,3] == "19:45411941:T:C",]
positions[positions[,3] == "19:45412079:C:T",]

# Get the APOE4 status
gt_rs429358 <- vcf@gt[positions[,3] == "19:45411941:T:C",]
gt_rs7412 <- vcf@gt[positions[,3] == "19:45412079:C:T",]

# Get info about the SNPs
info_rs429358 <- vcf@fix[positions[,3] == "19:45411941:T:C",]
info_rs7412 <- vcf@fix[positions[,3] == "19:45412079:C:T",]

# Save results
save(gt_rs429358, gt_rs7412, info_rs429358, info_rs7412, meta, file = "gt_results.RData")
