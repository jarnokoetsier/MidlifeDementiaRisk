
# Load packages
library(tidyverse)
library(ggpubr)

# Set working directory
setwd("C:/Users/Gebruiker/Documents/GitHub/Epi-LIBRA")

# Load meta data (methylation)
load("E:/Thesis/EXTEND/Phenotypes/metadat_extend.RData")
metaData <- dat

# Load genotype data
genotype_file <- read.table("E:/Thesis/EXTEND/Genotypes/RawData1/EXTEND_genotype_raw.fam")
genotype_file <- genotype_file[genotype_file$V2 %in% metaData$ID,]

# Load pre-processed methylation data
load("E:/Thesis/MLData/pheno.RData")

# Load processed genotype data
load("E:/Thesis/EXTEND/Genotypes/Plots/gt_results.RData")
gt_rs7412 <- gt_rs7412[names(gt_rs7412) %in% metaData$ID]

###############################################################################

# Plot age distribution

###############################################################################

# Prepare meta data object
metaData$AgeGroup <- rep("Mid-life",nrow(metaData))
metaData$AgeGroup[metaData$Age < 40] <- "Early-life"
metaData$AgeGroup[metaData$Age > 75] <- "Late-life"
metaData$AgeGroup <- factor(metaData$AgeGroup,
                            levels = c("Early-life",
                                       "Mid-life",
                                       "Late-life"))

# Make color palette
pal <- RColorBrewer::brewer.pal(n = 4, name = "Reds")[2:4]
names(pal) <- c("Early-life", "Mid-life", "Late-life")

x_min <- 40
x_max <- 75
# Make top of plot
top <- ggplot() +
  geom_rect(aes(xmin = min(metaData$Age)-5, xmax = x_min, ymin = 0, ymax = 1), fill = pal[1], alpha = 0.8) +
  geom_rect(aes(xmin = x_min, xmax = x_max, ymin = 0, ymax = 1), fill = pal[2], alpha = 0.8) +
  geom_rect(aes(xmin = x_max, xmax = max(metaData$Age) + 5, ymin = 0, ymax = 1), fill = pal[3], alpha = 0.8) +
  geom_vline(xintercept = x_min, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = x_max, linetype = "dashed", color = "black", size = 1) +
  geom_text(aes(x = min(metaData$Age) - 5 + (x_min - (min(metaData$Age) - 5))/2, y = 0.5, label = "Early-life"))+
  geom_text(aes(x = x_min + (x_max - x_min)/2, y = 0.5, label = "Mid-life"))+
  geom_text(aes(x = x_max + (max(metaData$Age) + 5 - x_max)/2, y = 0.5, label = "Late-life"))+
  theme_void()

# Make main plot
main <- ggplot() +
  geom_histogram(data = metaData, 
                 aes(x = Age, fill = AgeGroup), binwidth = 1, color = "grey") +
  #facet_grid(.~AgeGroup, scales = "free", space = "free") +
  geom_vline(xintercept = x_min, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = x_max, linetype = "dashed", color = "black", size = 1) +
  xlab("Age (years)") +
  ylab("Count") +
  xlim(c(min(metaData$Age) -5,max(metaData$Age) + 5)) +
  #ggtitle("Age Distribution") +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic"))

p <- ggarrange(top, main, nrow = 2, ncol = 1, heights = c(1,10), align = "v")

ggsave(p, file = "AgeDistribution.png", width = 8, height = 6)


###############################################################################

# Plot data overlap

###############################################################################

ageFiltered <- metaData[(metaData$Age <=75) & (metaData$Age >= 40),]
plotOverlap <- data.frame(ID = c(unique(metaData$ID),
                                 unique(genotype_file$V2),
                                 unique(pheno$IID),
                                 unique(ageFiltered$ID),
                                 unique(names(gt_rs7412)[-1])),
                          Source = c(rep("Methylation",length(unique(metaData$ID))),
                                     rep("Genotype", length(unique(genotype_file$V2))),
                                     rep("Methylation Processed",length(unique(pheno$IID))),
                                     rep("Methylation Age Filtered",length(unique(ageFiltered$ID))),
                                     rep("Genotype Processed",length(unique(names(gt_rs7412)[-1])))
                                     ))

orderID <- names(rev(sort(table(plotOverlap$ID[!str_detect(plotOverlap$Source,"Genotype")]))))
orderID <- names(rev(sort(table(plotOverlap$ID))))
plotOverlap$ID <- factor(plotOverlap$ID, levels = orderID)
plotOverlap$Source <- factor(plotOverlap$Source,
                             c("Methylation", "Methylation Age Filtered", "Methylation Processed",
                               "Genotype", "Genotype Processed"))


p <- ggplot() +
  geom_tile(data = plotOverlap, aes(x = ID, y = Source, fill = Source)) +
  geom_text(aes(x = orderID[round(length(orderID)/2)],
            y =  c("Methylation", "Methylation Age Filtered", "Methylation Processed",
                   "Genotype", "Genotype Processed"),
            label = c(length(unique(metaData$ID)), length(unique(ageFiltered$ID)),
                      length(unique(pheno$IID)), length(unique(genotype_file$V2)),
                      length(names(gt_rs7412)[-1]))),
            color = "white"
            ) +
  scale_fill_brewer(palette = "Set1") +
  xlab("Samples") +
  ylab(NULL) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_blank())

ggsave(p, file = "OverlapData.png", width = 8, height = 6)



