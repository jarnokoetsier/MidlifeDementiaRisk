# Load packages
library(readxl)
library(tidyverse)
library(ggpubr)

# Set working directory
setwd("E:/Thesis/EXTEND/Phenotypes")

# Load data
load("CAIDE.Rdata")
load("EPILIBRA.Rdata")

scoreAll <- inner_join(CAIDE[,c(1,10:16)], EPILIBRA[,c(1,10:20)], by = c("ID" = "ID"))
colnames(scoreAll) <- c("ID", "Age", "Sex", "Education", "Systolic Blood Pressure",
                        "BMI", "Total Serum Cholesterol", "Physical Inactivity",
                        "Diet", " Physical Inactivity", "Smoking", "Alcohol Intake",
                        "Obesity", "Depression", "Type 2 Diabetes", "Hypertension", "HDL Cholesterol", 
                        "Heart Disease", "Kidney Disease")
rownames(scoreAll) <- scoreAll$ID
scoreAll <- scoreAll[,-1]


test <- as.data.frame(corrr::correlate(scoreAll, diagonal = 1))
rownames(test) <- test$term
test <- test[,-1]

plotDF <- gather(test)
plotDF$Source <- rep(rownames(test), ncol(test))
plotDF$Score_key <- rep(c(rep("CAIDE1",7), rep("LIBRA",11)), each = ncol(test))
plotDF$Score_source <- rep(c(rep("CAIDE1",7), rep("LIBRA",11)), ncol(test))


# CAIDE vs CAIDE
caide_caide <- ggplot(plotDF[(plotDF$Score_source == "CAIDE1") & (plotDF$Score_key == "CAIDE1"),]) +
  geom_point(aes(x = key, y = Source, color = value, size = abs(value))) +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))
#ggsave(caide_caide, file = "caide_caide.png", height = 7, width = 7)

# LIBRA vs LIBRA
libra_libra <- ggplot(plotDF[(plotDF$Score_source == "LIBRA") & (plotDF$Score_key == "LIBRA"),]) +
  geom_point(aes(x = key, y = Source, color = value, size = abs(value))) +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))
#ggsave(libra_libra, file = "libra_libra.png", height = 11, width = 11)

# LIBRA vs CAIDE
libra_caide <- ggplot(plotDF[(plotDF$Score_source == "LIBRA") & (plotDF$Score_key == "CAIDE1"),]) +
  geom_point(aes(y = key, x = Source, color = value, size = abs(value))) +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))

# CAIDE vs LIBRA
caide_libra <- ggplot(plotDF[(plotDF$Score_source == "LIBRA") & (plotDF$Score_key == "CAIDE1"),]) +
  geom_point(aes(x = key, y = Source, color = value, size = abs(value))) +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))


p <- ggarrange(libra_caide, caide_caide, libra_libra, caide_libra,
          nrow = 2, ncol = 2, align = "hv")


ggsave(p, file = "ScoreCor.png", height = 11, width = 10)

