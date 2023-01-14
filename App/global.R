#=============================================================================#
# File: global.R
# Date: October 27, 2022										                                      
# Author: Jarno Koetsier                                                      
# Data: 'trainingData_filtered.RData', 'trainingData.RData', trainingClass.Rdata'
# 'testData.RData', "finalModel.RData', 'sampleInfo_filtered.RData', and
# 'featureInfo.RData'
#
# R version: 4.2.1 (getRversion())
# RStudio version: 2022.7.1.544 (RStudio.Version())
#=============================================================================#

# DISCLAIMER: The code has been developed using R version 4.2.1. Although the 
# code will likely work for other R versions as well, its functionality for 
# other R versions can not be guaranteed. 

# Clear working environment
rm(list = ls())

# Required packages
CRANpackages <- c("tidyverse",         # Data formatting and plotting
                  "shiny",             # Make App
                  "shinyWidgets",      # Widgets for app
                  "shinycssloaders",   # Loading figure
                  "shinythemes",       # Layout theme for app
                  "caret",             # Machine learning workflow
                  "glmnet",            # Elastic net
                  "plotly")            # Interactive plots


#Install (if not yet installed) and load the required packages: 
for (pkg in CRANpackages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, ask = FALSE)
  require(as.character(pkg), character.only = TRUE)
}

#library(shinydashboard)
#library(shinydashboardPlus)
#library(shinyBS)

# Set working directory
wd <- getwd()
