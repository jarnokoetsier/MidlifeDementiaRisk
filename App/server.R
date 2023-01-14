#=============================================================================#
# File: server.R
# Date: October 27, 2022										                                      
# Author: Jarno Koetsier                                                      
# Data: NA
#
# R version: 4.2.1 (getRversion())
# RStudio version: 2022.7.1.544 (RStudio.Version())
#=============================================================================#

# DISCLAIMER: The code has been developed using R version 4.2.1. Although the 
# code will likely work for other R versions as well, its functionality for 
# other R versions can not be guaranteed. 

#Server
server <- function(input, output, session) {
  
  #******************************************************************************#
  # 1) Preparation
  #******************************************************************************#
  
  # Welcome message
  sendSweetAlert(
    session = session,
    title = paste0("Hi ", Sys.getenv("USERNAME")),
    text = tags$span("Welcome to the ", tags$b("Class Prediction App"), "! ", tags$br(),
                     "If you need information about how to use the app, please click
                     on the ", tags$em("Information"), " tab."),
    type = "info",
    html = TRUE,
    btn_labels = "Yeah let's go!",
    showCloseButton = TRUE)
  
}