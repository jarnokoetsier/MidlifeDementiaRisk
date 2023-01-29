
# Requared packages
library(tidyverse)
library(DT)
library(shinyWidgets)
library(shinycssloaders)
library(shinyFiles)
library("PRSMultiTrait")
#PRSMultiTrait::installDependenciesAndData()

# Get manifest
getManifest()
Traits <- Manifest_env$Ref_gwas_manifest$short

