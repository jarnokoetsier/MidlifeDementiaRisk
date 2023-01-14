#=============================================================================#
# File: ui.R
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

# User interface
ui <- fluidPage(theme = shinytheme("spacelab"),
                setBackgroundColor("#343434"),
                tags$head(tags$style(HTML("
                           .navbar-nav {
                           float: none !important;
                           }
                           .navbar-nav > li:nth-child(5) {
                           float: right;
                           }
                           .my_style_1{ 
                             background-image: url(Background1.jpg);
                           }
                           
                           .my_style_1 { margin-top: -20px; }
                           
                           .my_style_1 { width: 100%; }
                           
                           .container-fluid { padding-left: 0; padding-right: 0; }
                           
                           .my_style_1 { position: absolute; left: 0; }
                           
                           "))),
              
  
  fluidPage(
    # Allow warning/information messages
    useSweetAlert(),
    
    navbarPage(title = "EPI-DEM", id = "navbar",
               
               tabPanel("Home", 
                        value = "home", 
                        icon = icon("fas fa-home"),
                        class = "my_style_1",
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        
                        fluidRow(
                          column(4, offset = 4, 
                                 align = "center", 
                                 style = "background-color:#FFFFFF;",
                                 
                                 br(),
                                 
                                 h1(strong(span(style = "color:#000000", "Welcome to EPI-DEM!"))),
                                 
                                 h5(span(style = "color:#000000", "Get started by uploading your data.")),
                                 
                                 br(),
                                 
                                 prettyRadioButtons(
                                   inputId = "selectInput",
                                   label = NULL,
                                   choices = c("Raw Data",
                                               "Normalized Data"),
                                   status = "danger",
                                   inline = TRUE
                                 ),
                                 
                                 #Enter database accession  
                                 fileInput(inputId = "uploadData", 
                                           label = NULL),
                                 
                                 #Use example data
                                 actionBttn(inputId = "example", 
                                            label = "Example",
                                            style = "jelly",
                                            color = "royal",
                                            icon("fas fa-mouse-pointer")),

                                 #Start the analysis
                                 actionBttn(inputId = "startAnalysis",
                                            label = "Predict!",
                                            style = "jelly",
                                            color = "primary",
                                            icon = icon("arrow-right")),
                                 
                                 br(),
                                 br(),
                                 br()
                          ) #Column
                        ), #FluidRow
                        
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br()),
               
               ################################################################
               
               # CAIDE1
               
               ################################################################
               
               tabPanel("CAIDE1", 
                        value = "caide1",
                        
                        tabsetPanel(id = "tabs_caide1",
                                    
                                    #***************************************************#
                                    # Score
                                    #***************************************************#
                                    
                                    tabPanel("Score", value = "caide1_score",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               
                                             ),
                                             
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                             )
                                             
                                    ),
                                    
                                    #***************************************************#
                                    # Category
                                    #***************************************************#
                                    
                                    tabPanel("Category", value = "caide1_category",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               
                                             ),
                                             
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                             )
                                    ),
                                    
                                    
                                    #***************************************************#
                                    # Factors
                                    #***************************************************#
                                    
                                    tabPanel("Factors" , value = "caide1_factors",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               
                                             )#End of mainPanel
                                             
                                    )
                                    
                        ) # End of tabset panel
                        
               ),#tabPanel
               
               
               ################################################################
               
               # CAIDE2
               
               ################################################################
               tabPanel("CAIDE2", 
                        value = "caide2",
                        
                        tabsetPanel(id="tabs_caide2",
                                    
                                    #***************************************************#
                                    # Score
                                    #***************************************************#
                                    
                                    tabPanel("Score", value = "caide2_score",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               
                                             ),
                                             
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                               
                                             )#mainPanel
                                             
                                    ),#tabPanel
                                    
                                    #***************************************************#
                                    # Category
                                    #***************************************************#
                                    tabPanel("Category",value = "caide2_category",
                                             
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                               
                                             ),#sidebarPanel
                                             mainPanel(
                                               
                                             )#mainPanel
                                             
                                    ),#tabPanel
                                    
                                    #***************************************************#
                                    # Factors
                                    #***************************************************#
                                    
                                    tabPanel("Factors" , value = "caide2_factors",
                                             br(),
                                             #==========================================#
                                             # Side bar
                                             #==========================================#
                                             sidebarPanel(
                                             ),
                                             
                                             #==========================================#
                                             # main Panel
                                             #==========================================#
                                             mainPanel(
                                               
                                             )#End of mainPanel
                                             
                                    )#tabPanel
                                    
                                    
                        )#tabSetPanel
               ),#tabPanel
               
               ################################################################
               #
               # LIBRA
               # 
               # ################################################################
               tabPanel("LIBRA",
                        value = "libra",
                        tabsetPanel(id="tabs_libra",
                                    
                                    #***************************************************#
                                    # Score
                                    #***************************************************#
                                    tabPanel("Score",value = "libra_score",
                                             br(),
                                             #==========================================#
                                             # Side bar panel
                                             #==========================================#
                                             sidebarPanel(
                                             ),
                                             
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                             )
                                    ),
                                    
                                    #***************************************************#
                                    # Category
                                    #***************************************************#
                                    tabPanel("Category", value = "libra_category",
                                             br(),
                                             #==========================================#
                                             # Side bar panel
                                             #==========================================#
                                             sidebarPanel(
                                             ),
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                             )
                                             
                                    ),
                                    #***************************************************#
                                    # Factors
                                    #***************************************************#
                                    tabPanel("Factors", value = "libra_factors",
                                             br(),
                                             #==========================================#
                                             # Side bar panel
                                             #==========================================#
                                             sidebarPanel(
                                             ),
                                             #==========================================#
                                             # Main panel
                                             #==========================================#
                                             mainPanel(
                                             )
                                             
                                    )
                                    
                                    
                        )#tabsetPanel
               ),#tabPanel
               
               ################################################################
               #
               # # Documentation
               # 
               # ################################################################
               tabPanel("Documentation",
                        value = "Documentation",
                        icon = icon("far fa-question-circle"),
                        tabsetPanel(id="tabs_documentation",
                                    
                                    #***************************************************#
                                    # Documentation of CAIDE1
                                    #***************************************************#
                                    tabPanel("CAIDE1",value = "doc_caide1"
                                    ),
                                    
                                    #***************************************************#
                                    # Documentation of CAIDE2
                                    #***************************************************#
                                    tabPanel("CAIDE2",value = "doc_caide2"
                                    ),
                                    
                                    #***************************************************#
                                    # Documentation of LIBRA
                                    #***************************************************#
                                    tabPanel("LIBRA",value = "doc_libra"
                                    ),
                        )# tabsetpanel
               ) #tabPanel
               
               
               
    ) #navbar page
    
  )#fluidPage
)#tagList