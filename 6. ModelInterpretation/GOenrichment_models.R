# Load packages
library(wateRmelon)
library(patchwork)
library(ggpubr)
library(rrvgo)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("~/allModels.RData")
load("~/Data/X_test.RData")

# How many features per model?
nFeatures <- rep(NA, length(allModels))
for (i in 1:length(allModels)){
  varImportance <- varImp(allModels[[i]])$importance
  selCpGs <- rownames(varImportance)[varImportance[,1] != 0]
  
  nFeatures[i] <- length(selCpGs)
}

age_cpgs <- read.csv("aging-10-101508-s005.csv")
selCpGs <- intersect(age_cpgs$ID[-1],rownames(X_test))
nFeatures <- c(nFeatures, length(selCpGs))

# Function to capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Perform GO overrepresentation analysis on the models' features
pvalues <- matrix(NA, nrow = 22708, ncol = length(allModels))
FDRs <- pvalues
for (i in 1:length(allModels)){
  varImportance <- varImp(allModels[[i]])$importance
  selCpGs <- rownames(varImportance)[varImportance[,1] != 0]
  allCpGs <- rownames(X_test)
  
  gst <- gometh(sig.cpg=selCpGs, all.cpg=allCpGs, collection="GO",
                array.type = "EPIC",
                plot.bias=TRUE)
  
  pvalues[,i] <- gst$P.DE
  FDRs[,i] <- gst$FDR
  
}
pvalues <- as.data.frame(pvalues)
FDRs <- as.data.frame(FDRs)
colnames(pvalues) <-  names(allModels)
rownames(pvalues) <- firstup(paste0(gst$TERM, " (", gst$ONTOLOGY, ")"))
colnames(FDRs) <- names(allModels)
rownames(FDRs) <- firstup(paste0(gst$TERM, " (", gst$ONTOLOGY, ")"))
FDRs <- as.data.frame(FDRs)

# Also for age: available at https://doi.org/10.18632%2Faging.101508
age_cpgs <- read.csv("aging-10-101508-s005.csv")
selCpGs <- intersect(age_cpgs$ID[-1],rownames(X_test))
allCpGs <- rownames(X_test)
gst <- gometh(sig.cpg=selCpGs, all.cpg=allCpGs, collection="GO",
              array.type = "EPIC",
              plot.bias=TRUE)
pvalues$Age<- gst$P.DE
FDRs$Age <- gst$FDR

# Which terms reach significance in more than 3 models?
test <- pvalues[rowSums(pvalues < 0.05) > 3,]
terms <- gst[rownames(test),]
colnames(test) <- c("BMI", "Type II Diabetes", "L-M Alcohol Intake",
                    "HDL Cholesterol", "Total Cholesterol", "Physical Inactivity",
                    "Heart Disease", "Education", "Depression", "Systolic Blood Pressure",
                    "Dietary Intake", "Sex", "Smoking", "Age")

# Prepare data for plotting
plotDF <- gather(test)
plotDF$Name <- rep(rownames(test), ncol(test))
plotDF$Sig <- ifelse(plotDF$value < 0.05, "Yes", "No")

# Perform hierarchical clustering on the terms
clusters <- hclust(dist(-log(test)), method = "ward.D2")
order <- clusters$labels[clusters$order]
plotDF$Name <- factor(plotDF$Name, levels = order)

# Make plot (heatmap)
p <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Name, fill = -log(value), color = Sig),
            width = 0.9, height = 0.9, linewidth = 0.5) +
  scale_fill_viridis_c() +
  scale_color_manual(values = c("white","black")) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = "-log p-value") +
  guides(color = "none") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

# Save plot
ggsave(p, file = "heatmap_GOenrichment_models.png", width = 10, height = 8)

# Save GO results
save(test, file = "GOenrichment_models.RData")


# Now reduce similar GO terms
load("GOenrichment_models.RData")

# Prepare data
gst$Name <- firstup(paste0(gst$TERM, " (", gst$ONTOLOGY, ")"))
gst$ID <- rownames(gst)
rownames(gst) <- gst$Name
gst_fil <- gst[rownames(test),]
test_copy <- test
rownames(test_copy) <- gst_fil$ID

# BP
simMatrix_BP <- calculateSimMatrix(gst_fil$ID[gst_fil$ONTOLOGY == "BP"],
                                orgdb = "org.Hs.eg.db",
                                ont = "BP", 
                                method = "Rel")

reduceTerms_BP <- reduceSimMatrix(simMatrix_BP,
                               rowMeans(-log(test_copy[gst_fil$ONTOLOGY == "BP",])),
                               threshold = 0.7,
                               orgdb = "org.Hs.eg.db")

# MF
simMatrix_MF <- calculateSimMatrix(gst_fil$ID[gst_fil$ONTOLOGY == "MF"],
                                   orgdb = "org.Hs.eg.db",
                                   ont = "MF", 
                                   method = "Rel")

reduceTerms_MF <- reduceSimMatrix(simMatrix_MF,
                                  rowMeans(-log(test_copy[gst_fil$ONTOLOGY == "MF",])),
                                  threshold = 0.7,
                                  orgdb = "org.Hs.eg.db")

# CC
simMatrix_CC <- calculateSimMatrix(gst_fil$ID[gst_fil$ONTOLOGY == "CC"],
                                   orgdb = "org.Hs.eg.db",
                                   ont = "CC", 
                                   method = "Rel")

reduceTerms_CC <- reduceSimMatrix(simMatrix_CC,
                                  rowMeans(-log(test_copy[gst_fil$ONTOLOGY == "CC",])),
                                  threshold = 0.7,
                                  orgdb = "org.Hs.eg.db")

# Combine reduced terms
reduceTerms <- rbind.data.frame(reduceTerms_BP, reduceTerms_MF, reduceTerms_CC)
reduceTerms <- reduceTerms[gst_fil$ID,]
all(gst_fil$Name == rownames(test))
all(gst_fil$ID == rownames(reduceTerms))
reduceTerms$Name<- gst_fil$Name
terms <- gst[rownames(test),]
colnames(test) <- c("BMI", "Type II Diabetes", "L-M Alcohol Intake",
                    "HDL Cholesterol", "Total Cholesterol", "Physical Inactivity",
                    "Heart Disease", "Education", "Depression", "Systolic Blood Pressure",
                    "Dietary Intake", "Sex", "Smoking", "Age")

# Prepare data for plotting
plotDF <- gather(test)
plotDF$Name <- rep(rownames(test), ncol(test))
plotDF$Sig <- ifelse(plotDF$value < 0.05, "Yes", "No")

# Perform hierarchical clustering
clusters <- hclust(dist(-log(test)), method = "ward.D2")
order <- clusters$labels[clusters$order]
plotDF$Name <- factor(plotDF$Name, levels = order)

# Main plot: heatmap
main <- ggplot(plotDF) +
  geom_tile(aes(x = key, y = Name, fill = -log(value), color = Sig),
            width = 0.9, height = 0.9, linewidth = 0.5) +
  scale_fill_viridis_c() +
  scale_color_manual(values = c("white","black")) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = "-log p-value") +
  guides(color = "none") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

# Side plot: parent terms
reduceTerms$Name <- factor(reduceTerms$Name, levels = order)

colors <- c(RColorBrewer::brewer.pal(n = 8, name = "Reds")[3:5], 
            RColorBrewer::brewer.pal(n = 8, name = "Blues")[3:6],
            RColorBrewer::brewer.pal(n = 8, name = "Purples")[3:6])

colors <- RColorBrewer::brewer.pal(n = 9, "Set3")
side <- ggplot(reduceTerms) +
  geom_tile(aes(x = 1, y = Name, fill = parent),
            width = 0.9, height = 0.9, linewidth = 0.5) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = c("white","black")) +
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

# Combine plots
p <- ggarrange(side, main, ncol = 2, nrow = 1, align = "h", widths = c(5.5,5))

# Save figure
ggsave(p, file = "heatmap_GOenrichment_models_groups.png", width = 10, height = 8)
