# Load packages
library(readxl)
library(tidyverse)
library(ggpubr)
library(patchwork)

# Set working directory
setwd("E:/Thesis/EXTEND/Phenotypes")

# Load data
load("CAIDE.Rdata")
load("CAIDE2.RData")
load("EPILIBRA.Rdata")

###############################################################################

# Correlation of CAIDE and LIBRA

###############################################################################

scoreAll <- inner_join(CAIDE2[,c(1,10:17)], EPILIBRA[,c(1,10:20)], by = c("ID" = "ID"))
colnames(scoreAll) <- c("ID", "Age", "Sex", "Education", "Systolic Blood Pressure",
                        "BMI", "Total Cholesterol", "Physical Inactivity", "APOE \u03b54 status",
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
plotDF$Score_key <- rep(c(rep("CAIDE",8), rep("LIBRA",11)), each = ncol(test))
plotDF$Score_source <- rep(c(rep("CAIDE",8), rep("LIBRA",11)), ncol(test))


# CAIDE vs CAIDE
caide_caide <- ggplot(plotDF[(plotDF$Score_source == "CAIDE") & (plotDF$Score_key == "CAIDE"),]) +
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
  ylab("LIBRA") +
  xlab("LIBRA") +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold",
                                    size = 12),
        axis.title.y = element_text(angle = 90, 
                                    vjust = 1,
                                    face = "bold",
                                    size = 12),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))
#ggsave(libra_libra, file = "libra_libra.png", height = 11, width = 11)

# LIBRA vs CAIDE
libra_caide <- ggplot(plotDF[(plotDF$Score_source == "LIBRA") & (plotDF$Score_key == "CAIDE"),]) +
  geom_point(aes(y = key, x = Source, color = value, size = abs(value))) +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  ylab("CAIDE") +
  xlab("LIBRA") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90, 
                                    vjust = 1,
                                    face = "bold",
                                    size = 12),
        legend.position = "none",
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))

# CAIDE vs LIBRA
caide_libra <- ggplot(plotDF[(plotDF$Score_source == "LIBRA") & (plotDF$Score_key == "CAIDE"),]) +
  geom_point(aes(x = key, y = Source, color = value, size = abs(value))) +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  xlab("CAIDE") +
  ylab("LIBRA") +
  theme_minimal() +
  theme(axis.title.x = element_text(face = "bold",
                                    size = 12),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 14),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10))



p <- libra_caide + caide_caide + libra_libra + caide_libra + 
  plot_layout(nrow = 2)



ggsave(p, file = "ScoreCor.png", height = 10, width = 10)

library(grid)
legendPlot <- ggplot(plotDF[(plotDF$Score_source == "LIBRA") & (plotDF$Score_key == "LIBRA"),]) +
  geom_point(aes(x = key, y = Source, color = value, size = abs(value))) +
  scale_color_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  scale_size_continuous(limits = c(0,1)) +
  ylab("LIBRA") +
  xlab("LIBRA") +
  labs(color = "Spearman\nCorrelation", size = "|Spearman\nCorrelation|") +
  theme_minimal()
legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)

###############################################################################

# Constitution of each score

###############################################################################

#*****************************************************************************#
# CAIDE1
#*****************************************************************************#

Value <- c(table(CAIDE$age_c)/nrow(CAIDE),
           table(CAIDE$Sex_c)/nrow(CAIDE),
           table(CAIDE$Edu_c)/nrow(CAIDE),
           table(CAIDE$Syst_c)/nrow(CAIDE),
           table(CAIDE$BMI_c)/nrow(CAIDE),
           table(CAIDE$Chol_c)/nrow(CAIDE),
           table(CAIDE$PHYSICAL_c)/nrow(CAIDE)
)

Score <- c(names(table(CAIDE$age_c)/nrow(CAIDE)),
           names(table(CAIDE$Sex_c)/nrow(CAIDE)),
           names(table(CAIDE$Edu_c)/nrow(CAIDE)),
           names(table(CAIDE$Syst_c)/nrow(CAIDE)),
           names(table(CAIDE$BMI_c)/nrow(CAIDE)),
           names(table(CAIDE$Chol_c)/nrow(CAIDE)),
           names(table(CAIDE$PHYSICAL_c)/nrow(CAIDE))
)

ScoreFactor <-  c(rep("Age",3), 
                  rep("Sex",2), 
                  rep("Education",2), 
                  rep("Systolic Blood Pressure",2),
                  rep("BMI",2), 
                  rep("Total Cholesterol",2), 
                  rep("Physical Inactivity",2)
)



plotDF <- data.frame(Value, Score, ScoreFactor)

colors <- c("grey",RColorBrewer::brewer.pal(n = 8, name = "Reds")[2:5])

p <- ggplot(plotDF) +
  geom_bar(aes(x = ScoreFactor, y = Value, fill = Score), 
           stat = "identity", color = 'black') +
  scale_fill_manual(values = colors) +
  xlab("") +
  ylab("Sample Proportion") +
  ggtitle("CAIDE1 Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

ggsave(p, file = "BarChart_Factors_CAIDE1.png", width = 8, height = 6)

#*****************************************************************************#
# CAIDE2
#*****************************************************************************#

Value <- c(table(CAIDE2$age_c)/nrow(CAIDE2),
           table(CAIDE2$Sex_c)/nrow(CAIDE2),
           table(CAIDE2$Edu_c)/nrow(CAIDE2),
           table(CAIDE2$Syst_c)/nrow(CAIDE2),
           table(CAIDE2$BMI_c)/nrow(CAIDE2),
           table(CAIDE2$Chol_c)/nrow(CAIDE2),
           table(CAIDE2$PHYSICAL_c)/nrow(CAIDE2),
           table(CAIDE2$APOE_c)/nrow(CAIDE2)
)

Score <- c(names(table(CAIDE2$age_c)/nrow(CAIDE2)),
           names(table(CAIDE2$Sex_c)/nrow(CAIDE2)),
           names(table(CAIDE2$Edu_c)/nrow(CAIDE2)),
           names(table(CAIDE2$Syst_c)/nrow(CAIDE2)),
           names(table(CAIDE2$BMI_c)/nrow(CAIDE2)),
           names(table(CAIDE2$Chol_c)/nrow(CAIDE2)),
           names(table(CAIDE2$PHYSICAL_c)/nrow(CAIDE2)),
           names(table(CAIDE2$APOE_c)/nrow(CAIDE2))
)

ScoreFactor <-  c(rep("Age",3), 
                  rep("Sex",2), 
                  rep("Education",2), 
                  rep("Systolic Blood Pressure",2),
                  rep("BMI",2), 
                  rep("Total Cholesterol",2), 
                  rep("Physical Inactivity",2),
                  rep("APOE \u03b54 status",2)
)



plotDF <- data.frame(Value, Score, ScoreFactor)

colors <- c("grey",RColorBrewer::brewer.pal(n = 8, name = "Reds")[c(2:4,6)])

p <- ggplot(plotDF) +
  geom_bar(aes(x = ScoreFactor, y = Value, fill = Score), 
           stat = "identity", color = 'black') +
  scale_fill_manual(values = colors) +
  xlab("") +
  ylab("Sample Proportion") +
  ggtitle("CAIDE2 Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

ggsave(p, file = "BarChart_Factors_CAIDE2.png", width = 8, height = 6)



#*****************************************************************************#
# LIBRA
#*****************************************************************************#

Value <- c(table(EPILIBRA$MEDITERANIAN)/nrow(EPILIBRA),
           table(EPILIBRA$PHYSICAL_INACTIVITY)/nrow(EPILIBRA),
           table(EPILIBRA$SMOKING)/nrow(EPILIBRA),
           table(EPILIBRA$LtoMAlcohol)/nrow(EPILIBRA),
           table(EPILIBRA$OBESITY)/nrow(EPILIBRA),
           table(EPILIBRA$DEPRESSION)/nrow(EPILIBRA),
           table(EPILIBRA$DIABETEII)/nrow(EPILIBRA),
           table(EPILIBRA$HYPERTENSTION)/nrow(EPILIBRA),
           table(EPILIBRA$Highcholesterol)/nrow(EPILIBRA),
           table(EPILIBRA$Heartdisease)/nrow(EPILIBRA),
           table(EPILIBRA$kidneydisease)/nrow(EPILIBRA)
)

Score <-  c(names(table(EPILIBRA$MEDITERANIAN)/nrow(EPILIBRA)),
            names(table(EPILIBRA$PHYSICAL_INACTIVITY)/nrow(EPILIBRA)),
            names(table(EPILIBRA$SMOKING)/nrow(EPILIBRA)),
            names(table(EPILIBRA$LtoMAlcohol)/nrow(EPILIBRA)),
            names(table(EPILIBRA$OBESITY)/nrow(EPILIBRA)),
            names(table(EPILIBRA$DEPRESSION)/nrow(EPILIBRA)),
            names(table(EPILIBRA$DIABETEII)/nrow(EPILIBRA)),
            names(table(EPILIBRA$HYPERTENSTION)/nrow(EPILIBRA)),
            names(table(EPILIBRA$Highcholesterol)/nrow(EPILIBRA)),
            names(table(EPILIBRA$Heartdisease)/nrow(EPILIBRA)),
            names(table(EPILIBRA$kidneydisease)/nrow(EPILIBRA))
)

ScoreFactor <-  c(rep("Diet",2), 
                  rep("Physical Inactivity",2), 
                  rep("Smoking",2), 
                  rep("Alcohol Intake",2),
                  rep("Obesity",2), 
                  rep("Depression",2), 
                  rep("Type 2 Diabetes",2),
                  rep("Hypertension",2),
                  rep("HDL Cholesterol",2),
                  rep("Heart Disease",2),
                  rep("Kidney Disease",2)
)



plotDF <- data.frame(Value, Score, ScoreFactor)
plotDF$Score <- factor(plotDF$Score,levels = as.character(sort(as.numeric(unique(plotDF$Score)))))

colors <- c(RColorBrewer::brewer.pal(n = 3, name = "Blues")[c(3,2)], "grey",
            RColorBrewer::brewer.pal(n = 8, name = "Reds")[2:8])
p <- ggplot(plotDF) +
  geom_bar(aes(x = ScoreFactor, y = Value, fill = Score), 
           stat = "identity", color = 'black') +
  scale_fill_manual(values = colors) +
  xlab("") +
  ylab("Sample Proportion") +
  ggtitle("LIBRA Score") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

ggsave(p, file = "BarChart_Factors_LIBRA.png", width = 8, height = 6)

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


###############################################################################

# Distribution of each score

###############################################################################

p <- ggplot(CAIDE) +
  geom_bar(aes(x = CAIDE, y = after_stat(count)/919), 
           fill = "#DC3535", color = "black", linewidth = 0.8) +
  xlab("CAIDE1 Score") +
  ylab("Sample Proportion")  +
  ggtitle("CAIDE1 Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))
ggsave(p, file = "Distribution_CAIDE1.png", width = 8, height = 6)

p <- ggplot(CAIDE2) +
  geom_bar(aes(x = CAIDE2, y = after_stat(count)/789), 
           fill = "#F49D1A", color = "black", linewidth = 0.8) +
  xlab("CAIDE2 Score") +
  ylab("Sample Proportion") +
  ggtitle("CAIDE2 Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))
ggsave(p, file = "Distribution_CAIDE2.png", width = 8, height = 6)


p <- ggplot(EPILIBRA) +
  geom_histogram(aes(x = LIBRA, y = after_stat(count)/883), 
                 bins = 12, fill = "#B01E68", color = "black", linewidth = 0.8) +
  xlab("LIBRA Score") +
  ylab("Sample Proportion") +
  ggtitle("LIBRA Score") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))
ggsave(p, file = "Distribution_LIBRA.png", width = 8, height = 6)



