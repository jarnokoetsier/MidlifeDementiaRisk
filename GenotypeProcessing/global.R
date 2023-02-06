
# Requared packages
library(tidyverse)
library(corrr)
library(ggdendroplot)
library(patchwork)
library(DT)
library(shinyWidgets)
library(shinycssloaders)
library(shinyFiles)
library("PRSMultiTrait")
#PRSMultiTrait::installDependenciesAndData()

# Get manifest
getManifest()
Traits <- Manifest_env$Ref_gwas_manifest$short


