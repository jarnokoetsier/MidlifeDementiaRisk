
# Load packages
library(vcfR)
library(tidyverse)
library(ggpubr)
library(grid)

# Set working directory
setwd("E:/Thesis/EXTEND/Genotypes/Plots")

################################################################################

# Read VCF file

################################################################################

# Read vcf file of chromosome 19
vcf <- read.vcfR(file = "chr19.dose.vcf.gz")

# Investigate structure of vcf object
head(vcf@gt[,1:5])
head(vcf@fix)

# Get positions of SNPs
positions <- vcf@fix

# Get meta data 
meta <- vcf@meta

# APOE status is determined by two SNPs:
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


################################################################################

# Extract APOE status

################################################################################

# Load data
load("gt_results.RData")

# Check if samples are in same order
all(names(gt_rs429358) == names(gt_rs7412))

# Extract genotypes
GT_rs429358 <- unlist(lapply(str_split(gt_rs429358, ":"), `[[`, 1))[-1]
GT_rs429358 <- str_replace_all(GT_rs429358, "0", "T")
GT_rs429358 <- str_replace_all(GT_rs429358, "1", "C")

GT_rs7412 <- unlist(lapply(str_split(gt_rs7412, ":"), `[[`, 1))[-1]
GT_rs7412 <- str_replace_all(GT_rs7412, "0", "C")
GT_rs7412 <- str_replace_all(GT_rs7412, "1", "T")


# Get haplotypes (TT, CC, TC, and CT)
hap1 <- paste0(unlist(lapply(str_split(GT_rs429358, "|"), `[[`, 2)),
               unlist(lapply(str_split(GT_rs7412, "|"), `[[`, 2)))

hap2 <- paste0(unlist(lapply(str_split(GT_rs429358, "|"), `[[`, 4)),
               unlist(lapply(str_split(GT_rs7412, "|"), `[[`, 4)))

# Use haplotypes to get the APOE status
e1 <- str_detect(hap1, "CT") + str_detect(hap2, "CT")
e2 <- str_detect(hap1, "TT") + str_detect(hap2, "TT")
e3 <- str_detect(hap1, "TC") + str_detect(hap2, "TC")
e4 <- str_detect(hap1, "CC") + str_detect(hap2, "CC")

# Collect into data frame
output <- data.frame(sampleID = names(gt_rs429358)[-1],
                     e1 = e1,
                     e2 = e2,
                     e3 = e3,
                     e4 = e4)

# Last check
all(rowSums(output[,2:5]) == 2)

# Write output file
write.csv(output, file = "APOE_status.csv", row.names = FALSE)


################################################################################

# Make plots

################################################################################

# Format for plotting
colnames(output) <- c("sampleID", "APOE \u03b51", "APOE \u03b52", "APOE \u03b53", "APOE \u03b54")
plotDF <- gather(output[,2:5])
plotDF$SampleID <- rep(output$sampleID,4)

e4_samples <- c(output$sampleID[output$`APOE ε4` == 2], 
                output$sampleID[(output$`APOE ε4` == 1) & (output$`APOE ε3` == 1)],
                output$sampleID[(output$`APOE ε4` == 1) & (output$`APOE ε2` == 1)],
                output$sampleID[(output$`APOE ε4` == 1) & (output$`APOE ε1` == 1)],
                output$sampleID[output$`APOE ε3` == 2],
                output$sampleID[(output$`APOE ε3` == 1) & (output$`APOE ε2` == 1)],
                output$sampleID[(output$`APOE ε3` == 1) & (output$`APOE ε1` == 1)],
                output$sampleID[output$`APOE ε2` == 2],
                output$sampleID[(output$`APOE ε2` == 1) & (output$`APOE ε1` == 1)],
                output$sampleID[output$`APOE ε1` == 2]
                )
plotDF$SampleID <- factor(plotDF$SampleID, levels = e4_samples)

# make plot
apoe_plot <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = SampleID, fill = as.factor(value))) +
  facet_grid(cols = vars(key), scales = "free", space = "free") +
  ylab("Samples") +
  #ggtitle("APOE status") +
  labs(fill = "# Alleles") +
  scale_fill_brewer(palette = "Reds") +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        strip.background = element_blank(),
        strip.text.x = element_blank())

# Save plot
ggsave(apoe_plot, file = "APOE_Status.png", width = 8, height = 6)
