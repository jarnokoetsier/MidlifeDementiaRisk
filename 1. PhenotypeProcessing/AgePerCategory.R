# Clear workspace and console
rm(list = ls())
cat("\014") 
gc()

# Load packages
library(readxl)
library(tidyverse)
library(ggpubr)
library(patchwork)

# Set working directory
setwd("E:/Thesis/EXTEND/Phenotypes")

# Load CAIDE data
load("CAIDE.Rdata")

# Divide CAIDE into three categories
CAIDE$Cat <- rep("Intermediate Risk", nrow(CAIDE))
CAIDE$Cat[CAIDE$CAIDE < 4] <- "Low Risk"
CAIDE$Cat[CAIDE$CAIDE > 7] <- "High Risk"
CAIDE$Cat <- factor(CAIDE$Cat, levels = c("Low Risk",
                                          "Intermediate Risk",
                                          "High Risk"))

# Plot the age per CAIDE category
p <- ggplot(CAIDE) +
  geom_rect(ymin = -Inf, ymax = 47, xmin = -Inf, xmax = Inf,fill = "#FFFFFF", alpha = 0.5) +
  geom_rect(ymin = 47, ymax = 53, xmin = -Inf, xmax = Inf,fill = "#F0F0F0", alpha = 0.5) +
  geom_rect(ymin = 53, ymax = Inf, xmin = -Inf, xmax = Inf,fill = "#D9D9D9" , alpha = 0.5) +
  geom_hline(yintercept = 47, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 53, linewidth = 1.5, linetype = "dashed", color = "black") +
  geom_boxplot(aes(x = Cat, y = Age, fill = Cat), 
               alpha = 1,outlier.shape = NA) +
  geom_point(aes(x= Cat, y = Age, color = Cat), size = 1.5, alpha = 0.8,
             position=position_jitterdodge(jitter.width = 0.7, jitter.height = 0)) +
  xlab("CAIDE1 Category") +
  ylab("Age") +
  scale_fill_manual(values = c("#FEE0D2","#FCBBA1", "#FC9272")) +
  scale_color_manual(values = c("#FB6A4A","#EF3B2C", "#CB181D")) +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic")) 

# Save plot
ggsave(p, file = "CAIDE1_Boxplot_Age.png", height = 3.9, width = 7)
