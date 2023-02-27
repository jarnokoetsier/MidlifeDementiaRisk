# Clear workspace and console
rm(list = ls())
cat("\014") 

# Install package
devtools::install_github("Rrtk2/PRS-multi-trait/Package/PRSMultiTrait")
library("PRSMultiTrait")
PRSMultiTrait::installDependenciesAndData()

library("PRSMultiTrait")
library(data.table)

files <- c("CKD_transEthnic_2019_Wuttke_curated.summaries",
           "CKD_European_2019_Wuttke_curated.summaries",
           "EA_2022_Okbay_curated.summaries",
           "CAD_2017_Nelson_curated.summaries",
           "PD_2019_Nalls_curated.summaries",
           "DC2_2020_Niarchou_curated.summaries")

short <- c("CKD_transEthnic",
           "CKD_European",
           "EA22",
           "CAD",
           "PD",
           "DC2")

year <- c("2019",
           "2019",
           "2022",
           "2017",
           "2019",
           "2020")

doi <- c("https://doi.org/10.1038%2Fs41588-019-0407-x",
         "https://doi.org/10.1038%2Fs41588-019-0407-x",
         "https://doi.org/10.1038/s41588-022-01016-z",
         "https://doi.org/10.1038/ng.3913", 
         "https://doi.org/10.1016/s1474-4422(19)30320-5",
         "https://www.nature.com/articles/s41467-018-04951-w",
         "https://www.nature.com/articles/s41398-020-0688-y#Sec2")


trait <- c("CKD_transEthnic (chronic kidney disease)",
           "CKD_European (chronic kidney disease)",
           "EA (Educational attainment)",
           "CAD (coronary artery disease)",
           "PD (Parkinson's disease)",
           "DC2 (Fish- and plant-based diet)")


type <- c("CAT",
          "CAT",
          "CONT",
          "CAT",
          "CAT",
          "CONT")
addGWAStoManifest()
for (i in 5:6){
  # load file
  dataObj <- fread(paste0("E:/Thesis/GWAS/New/", files[i]))
  sampleSize <- median(dataObj$N)
  
  # Add GWAS
  test(short = short[i], 
                                   n = sampleSize, 
                                   filename = paste0("E:/Thesis/GWAS/New/", files[i]), 
                                   year = year[i], 
                                   trait = trait[i], 
                                   DOI = doi[i], 
                                   genomeBuild = c("GRCh37"), 
                                   traitType = type[i], 
                                   rawSNPs = c("?"), 
                                   finalModelSNPs = c("?"), 
                                   modelRunningTime = c("?"),
                                   usedRefSet = c("?"), 
                                   processed = c(0), 
                                   FORCE = TRUE) 
  rm(dataObj)
  
  # prepare GWAS
  PRSMultiTrait::prepareGWAS(trait = short[i])
  
  # generation of PGM
  PRSMultiTrait::calcPGS_LDAK(Trait = short[i],Model = "bayesr")
}





test <- function (short = c("UniqueTraitName"), n = c(10000), filename = c("?"), 
          year = c("?"), trait = c("?"), DOI = c("?"), genomeBuild = c("?"), 
          traitType = c("?"), rawSNPs = c("?"), finalModelSNPs = c("?"), 
          modelRunningTime = c("?"), usedRefSet = c("?"), processed = c(0), 
          FORCE = FALSE) 
{
  getManifest()
  temp_man = data.frame(short = short, n = n, filename = filename, 
                        year = year, trait = trait, DOI = DOI, genomeBuild = genomeBuild, 
                        traitType = traitType, rawSNPs = rawSNPs, finalModelSNPs = finalModelSNPs, 
                        modelRunningTime = modelRunningTime, usedRefSet = usedRefSet, 
                        processed = processed)
  apply(temp_man, 2, function(x) {
    cat(paste0(x, "\n"))
  })
  if (!FORCE) {
    cat("\n>> Press [y] if information correct, then press [ENTER] <<")
    check = readline()
    if (check != "y") {
      return(message("Adding GWAS to Manifest aborted."))
    }
  }
  if (temp_man$short %in% Manifest_env$Ref_gwas_manifest$short) {
    message("Adding GWAS to Manifest failed!")
    return(message("'short' name is taken! Please check the input or make another unique name."))
  }
  Manifest_env$Ref_gwas_manifest[dim(Manifest_env$Ref_gwas_manifest)[1] + 
                                   1, ] = temp_man
  saveManifest()
  getManifest(1)
}

listofStandardizedGWASes <- short
