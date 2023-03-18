# Load package
library(wateRmelon)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("~/allModels.RData")
load("~/Data/X_test.RData")

# Function to capetilize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


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

# Add age: available at https://doi.org/10.18632%2Faging.101508
age_cpgs <- read.csv("aging-10-101508-s005.csv")
selCpGs <- intersect(age_cpgs$ID[-1],rownames(X_test))
allCpGs <- rownames(X_test)

gst <- gometh(sig.cpg=selCpGs, all.cpg=allCpGs, collection="GO",
              array.type = "EPIC",
              plot.bias=TRUE)

pvalues$Age<- gst$P.DE
FDRs$Age <- gst$FDR

test <- pvalues[rowSums(pvalues < 0.05) > 3,]

terms <- gst[rownames(test),]

colnames(test) <- c("BMI", "Type II Diabetes", "L-M Alcohol Intake",
                    "HDL Cholesterol", "Total Cholesterol", "Physical Inactivity",
                    "Heart Disease", "Education", "Depression", "Systolic Blood Pressure",
                    "Dietary Intake", "Sex", "Smoking", "Age")

plotDF <- gather(test)
plotDF$Name <- rep(rownames(test), ncol(test))
plotDF$Sig <- ifelse(plotDF$value < 0.05, "Yes", "No")


clusters <- hclust(dist(-log(test)), method = "ward.D2")
order <- clusters$labels[clusters$order]
plotDF$Name <- factor(plotDF$Name, levels = order)

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

ggsave(p, file = "heatmap_GOenrichment_models.png", width = 10, height = 8)
save(test, file = "GOenrichment_models.RData")
