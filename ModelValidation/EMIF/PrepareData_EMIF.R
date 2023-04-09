library(tidyverse)
library(caret)

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/metaData_EMIF.RData")
load("EMIF/predictedScore_factors_EMIF.RData")

#==============================================================================#
# Filter data
#==============================================================================#
metaData_EMIF <- as.data.frame(metaData_EMIF)
rownames(metaData_EMIF) <- metaData_EMIF$X

# Remove individuals with unmatching sex
metaData_EMIF <- metaData_EMIF[metaData_EMIF$sex.match == 1,]

# Keep midlife samples only
metaData_EMIF <- metaData_EMIF[metaData_EMIF$Age <= 75,]

# Remove converters
converters1 <-  unique(rownames(metaData_EMIF)[metaData_EMIF$CTR_Convert == 1])[-1]
converters2 <- unique(metaData_EMIF$X[(metaData_EMIF$LastFU_Diagnosis == "MCI") | (metaData_EMIF$LastFU_Diagnosis == "AD")])[-1]
converters2 <- intersect(converters2, metaData_EMIF$X[metaData_EMIF$Diagnosis == "NL"])
converters_all <- unique(c(converters1, converters2))
metaData_EMIF <- metaData_EMIF[setdiff(rownames(metaData_EMIF), converters_all),]
table(metaData_EMIF$CTR_Convert)

# Remove SCI individuals
SCI <- unique(rownames(metaData_EMIF)[metaData_EMIF$Diagnosis == "SCI"])
metaData_EMIF <- metaData_EMIF[setdiff(rownames(metaData_EMIF), SCI),]
table(metaData_EMIF$Diagnosis)

# Save meta data
samples <- intersect(metaData_EMIF$X, rownames(predictedScore_factors))
metaData_fil <- metaData_EMIF[samples,]
save(metaData_fil, file = "EMIF/metaData_fil.RData")

# save predicted scores
predictedScore_factors_fil <- predictedScore_factors[samples,]
save(predictedScore_factors_fil, file = "EMIF/predictedScore_factors_fil.RData")


table(metaData_fil$Diagnosis)








