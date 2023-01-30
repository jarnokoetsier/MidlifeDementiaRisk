
library(data.table)
bimFile <- fread("E:/Thesis/EXTEND/Genotypes/ChrBPData/Output_all/FINAL/EXTEND_PostImpute_FINAL_bp_dup.bim")

dupVar1 <- bimFile$V2[duplicated(bimFile$V2)]

setwd("E:/Thesis/EXTEND/Genotypes/")
write.table(dupVar, file = "plink.dupvar", quote = FALSE, col.names = FALSE, row.names = FALSE,
           sep = "\t")

devtools::install_github("Rrtk2/PRS-multi-trait/Package/PRSMultiTrait")




library("PRSMultiTrait")
#PRSMultiTrait::installDependenciesAndData()

getManifest()
Traits <- Manifest_env$Ref_gwas_manifest$short


for (t in Traits){
  predPRS(bfile = wslPath("E:/Thesis/EXTEND/Genotypes/ChrBPData/Output_all/FINAL/EXTEND_PostImpute_FINAL_bp_dup"), 
          Trait = t, 
          OverlapSNPsOnly=FALSE, 
          Force = FALSE)
}



famFile <- fread("E:/Thesis/EXTEND/Genotypes/ChrBPData/Output_all/FINAL/EXTEND_PostImpute_FINAL_bp_dup.fam")
df_Result_PGS = collect_all_PRS(cohort = "EXTEND_PostImpute_FINAL_bp_dup")
rownames(df_Result_PGS) <- famFile$V2
save(df_Result_PGS,file = "df_Result_PGS.RData")

library(corrr)
library(ggdendroplot)
library(patchwork)
# Calculate correlations
corrDF <- as.data.frame(correlate(df_Result_PGS, diagonal = 1, method = "spearman"))
rownames(corrDF) <- corrDF$term
corrDF <- corrDF[,-1]

# Format data for plotting
plotCor <- gather(corrDF)
plotCor$key1 <- rep(rownames(corrDF), ncol(corrDF))

# Perform clustering to get sample order
model <- hclust(as.dist(abs(1-corrDF)), "ward.D2")
order <- model$labels[model$order]
plotCor$key <- factor(plotCor$key, levels = order)
plotCor$key1 <- factor(plotCor$key1, levels = order)

main <- ggplot(plotCor) +
  geom_tile(aes(x = key, y = key1, fill = value)) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                        limits = c(-1,1)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = "Pearson\nCorrelation")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

dendroPlot <- ggplot() +
  geom_tile(data = plotCor, aes(x = as.numeric(key), y = 1), fill = "white", alpha = 0) +
  geom_dendro(model,xlim = c(1,length(unique(plotCor$key)))) +
  theme_void() +
  theme(legend.position = "none") 

dendroPlot + main +
  plot_layout(nrow = 2, ncol = 1,
              heights = c(1,3))

# Set working directory
setwd("E:/Thesis/EXTEND/Methylation")

# Load phenotype data
files <- list.files('Y')
for (f in files){
  load(paste0("Y/",f))
}

load("E:/Thesis/EXTEND/Phenotypes/metaData_ageFil.RData")

rownames(dat) <- dat$ID
total <- dat[rownames(df_Result_PGS),]

cor(df_Result_PGS$T2D, total$T2.Diabetes,use = "pairwise.complete.obs")
cor(df_Result_PGS$BMI, as.numeric(total$BMI),use = "pairwise.complete.obs")


trainingDF <- df_Result_PGS[Y_CAIDE2$ID,]

test <- cbind(trainingDF,Y_CAIDE2$CAIDE2)
colnames(test) <- c(colnames(test)[1:11], "Y")

model <- lm(data = test, Y ~.)

cor(test$ObesityClass1, Y_CAIDE2$BMI_c, method = "spearman")
