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

###############################################################################

# Constitution of each score

###############################################################################

library(tidyverse)

# Set working directory
setwd("E:/Thesis/EXTEND/Phenotypes")

# Load data
load("CAIDE.Rdata")
load("EPILIBRA.Rdata")

CAIDE1 <- CAIDE[,c(1,10:17)]
CAIDE1$Age3 <- ifelse(CAIDE1$age_c == 3, 3, 0)
CAIDE1$Age4 <- ifelse(CAIDE1$age_c == 4, 4, 0)
CAIDE1 <- CAIDE1[,c(1,2,10,11,3:9)]
colnames(CAIDE1) <- c("ID", "Age", "Age (>47)", "Age (>53)", "Sex", "Edu", "BP", "BMI", "Chol", "Physical", "CAIDE")
scores <- table(CAIDE1$CAIDE)

plotDF <- NULL
for (i in 1:length(scores)){

  temp <- CAIDE1[CAIDE1$CAIDE == names(scores)[i],]
  temp1 <- colSums(temp[,3:10] != 0)/scores[i]
  
  collect <- data.frame(Factor = names(temp1),
                        nonZero = temp1,
                        Score = rep(names(scores[i]), length(temp1)))
  
  plotDF <- rbind.data.frame(plotDF, collect)
}

plotDF$Score <- factor(plotDF$Score,
                       levels = as.character(0:12))

main <- ggplot(plotDF) +
  geom_bar(aes(y = nonZero, x = Score, fill = Factor), stat="identity") +
  facet_grid(rows = vars(Factor)) +
  xlab("CAIDE1 Score") +
  ylab("Proportion with non-zero score") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), labels = c("", "", "0.5", "", "1")) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text.x = element_text(),
        legend.position = "none",
        strip.background = element_rect(fill= "grey", linewidth = 0, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"))

CAIDE1$CAIDE <- factor(CAIDE1$CAIDE,
                       levels = as.character(0:12))
top <- ggplot(CAIDE1) +
  geom_bar(aes(x = CAIDE)) +
  ylab("# Samples") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

sidePlot <- data.frame(Factor = names( apply(CAIDE1[,3:10], 2, max)),
                       Value =  apply(CAIDE1[,3:10], 2, max))

side <- ggplot(sidePlot) +
  geom_bar(aes(y = Value, x = Factor), stat="identity") +
  ylab("Score Weight") +
  coord_flip() +
  facet_grid(rows = vars(Factor), scales = "free", space = "free") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(),
        #axis.line.y = element_blank(),
        #axis.ticks.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.y = element_blank())
  
  
  

library(patchwork)

p <- top + main + plot_spacer() + side +
  plot_layout(nrow = 2, byrow = FALSE) +
  plot_layout(heights = c(1,6), widths = c(6,1))

ggsave(p, file = "test.png", width = 8, height = 8)