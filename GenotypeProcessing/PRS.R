devtools::install_github("Rrtk2/PRS-multi-trait/Package/PRSMultiTrait")
library("PRSMultiTrait")
PRSMultiTrait::installDependenciesAndData()

getManifest(1)
predPRS(bfile = wslPath("E:/Thesis/EXTEND/Genotypes/ChrBPData/Output_all/FINAL/EXTEND_PostImpute_FINAL_bp_dup"), 
        Trait = "HDL", 
        OverlapSNPsOnly=FALSE, 
        Force = FALSE)