library(prospectr)
library(tidyverse)
library(caret)
library(patchwork)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/metaData_EMIF.RData")
load("EMIF/predictedScore_factors_EMIF.RData")

metaData_EMIF <- as.data.frame(metaData_EMIF)
rownames(metaData_EMIF) <- metaData_EMIF$X

# Remove individuals with unmatching sex
metaData_EMIF <- metaData_EMIF[metaData_EMIF$sex.match == 1,]

# Keep midlife samples only
metaData_EMIF <- metaData_EMIF[metaData_EMIF$Age <= 75,]

# Remove converters
converters1 <-  unique(rownames(metaData_EMIF)[metaData_EMIF$CTR_Convert == 1])[-1]
converters1 <- intersect(converters1, metaData_EMIF$X[(metaData_EMIF$Diagnosis == "NL")])
converters2 <- unique(metaData_EMIF$X[(metaData_EMIF$LastFU_Diagnosis == "MCI") | (metaData_EMIF$LastFU_Diagnosis == "AD")])[-1]
converters2 <- intersect(converters2, metaData_EMIF$X[(metaData_EMIF$Diagnosis == "NL")])
converters_all <- unique(c(converters1, converters2))
metaData_EMIF <- metaData_EMIF[setdiff(rownames(metaData_EMIF), converters_all),]
table(metaData_EMIF$CTR_Convert)
table(metaData_EMIF$Diagnosis)

samples <- intersect(rownames(predictedScore_factors), rownames(metaData_EMIF))
predictedScore_factors_fil <- predictedScore_factors[samples,]
metaData_fil <- metaData_EMIF[samples,]

all(metaData_fil$X == rownames(predictedScore_factors_fil))
table(metaData_fil$CTR_Convert)
table(metaData_fil$Diagnosis)

# Predict LIBRA and CAIDE
load("EMIF/Fit_CombineFactors_CAIDE1_RF.RData")
EpiCAIDE1 <- predict(fit, predictedScore_factors_fil)
load("EMIF/Fit_CombineFactors_LIBRA_RF.RData")
EpiLIBRA <- predict(fit, predictedScore_factors_fil)


plot_LIBRA <- data.frame(Diagnosis = metaData_fil$Diagnosis,
                         Value = EpiLIBRA,
                         Score = rep("Epi-LIBRA", length(EpiLIBRA)))

plot_CAIDE <- data.frame(Diagnosis = metaData_fil$Diagnosis,
                         Value = EpiCAIDE1,
                         Score = rep("Epi-CAIDE1", length(EpiCAIDE1)))

plot_all <- rbind.data.frame(plot_LIBRA, plot_CAIDE)
plot_all$Diagnosis[plot_all$Diagnosis == "NL"] <- "Control"
plot_all$Diagnosis <- factor(plot_all$Diagnosis, levels = c("Control", "SCI", "MCI", "AD"))

p <- ggplot(plot_all) +
  geom_boxplot(aes(x = Diagnosis, y = Value, fill = Score), alpha = 0.5) +
  geom_jitter(aes(x = Diagnosis, y = Value, color = Score), width = 0.15) +
  ylab("Risk Factor Score") +
  facet_grid(rows = vars(Score), scale = "free") +
  scale_fill_manual(values = c("#EF3B2C", "#807DBA")) +
  scale_color_manual(values = c("#EF3B2C", "#807DBA")) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(p, file = "Boxplot_LIBRA_CAIDE1.png", width = 6, height = 4)


test <- kruskal.test(Value ~ Diagnosis, data = plot_LIBRA)
test
pairwise.wilcox.test(plot_LIBRA$Value, plot_LIBRA$Diagnosis,
                     p.adjust.method = "BH")


test <- kruskal.test(Value ~ Diagnosis, data = plot_CAIDE)
test
pairwise.wilcox.test(plot_CAIDE$Value, plot_CAIDE$Diagnosis,
                     p.adjust.method = "BH")