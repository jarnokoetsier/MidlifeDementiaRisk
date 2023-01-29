
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
