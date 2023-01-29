ui <- tagList(
  
  tags$head(tags$style(HTML("
                           .navbar-nav {
                           float: none !important;
                           }
                           .navbar-nav > li:nth-child(6) {
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
    
    useSweetAlert(),
    
    navbarPage(title = "-   PRS Multi-Trait", id = "navbar",
               
               ###################################################################
               #  Data selection                                                   
               ###################################################################
               tabPanel("Data accession", 
                        value = "input_panel", 
                        icon = icon("fas fa-home"), class = "my_style_1",
                        
                        
                        #********************************************************#
                        #   Data selection
                        #********************************************************#
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
                                 
                                 h1(strong(span(style = "color:#000000", "Welcome to PRS Multi-Trait!"))),
                                 
                                 h5(span(style = "color:#000000", "Get started by uploading your PLINK files.")),
                                 
                                 br(),
                                 
                                 
                                 #Upload CELs
                                 shinyFilesButton(id = "bfile", 
                                                  label = "Click here to select a PLINK file",
                                                  title = "Please select a file:",
                                                  multiple = FALSE),
                                 br(),
                                 br()
                          )
                        ),
                        
                        
                        fluidRow(
                          
                          column(4, offset = 4, 
                                 align = "center", 
                                 style = "background-color:#FFFFFF;",
                                 
                                 #Use example data
                                 actionBttn(inputId = "example", 
                                            label = "Example",
                                            style = "jelly",
                                            color = "danger",
                                            icon("fas fa-mouse-pointer")),
                                 
                                 #Start the analysis
                                 actionBttn(inputId = "startAnalysis",
                                            label = "Start",
                                            style = "jelly",
                                            color = "danger",
                                            icon = icon("arrow-right")),
                                 
                                 
                                 
                                 br(),
                                 br(),
                                 
                          )
                        ),
                        
                        fluidRow(
                          column(4, offset = 4, align = "center",
                                 style = "background-color:#FFFFFF;",
                                 
                                 br(),
                                 awesomeCheckbox(inputId = "all_traits",
                                                 label = "Calculate PRS for all available traits",
                                                 value = TRUE,
                                                 status = "danger"),
                                 
                                 conditionalPanel(
                                   condition = "input.all_traits==false",
                                   
                                   selectInput(inputId = "traits_input",
                                               label = "Select traits",
                                               choices = Traits,
                                               multiple = TRUE),
                                 ),
                                 br()
                          )
                          
                        ),
                        
                        
                        #********************************************************#
                        #   Continue with saved data
                        #********************************************************#
                        
                        fluidRow(
                          column(4, offset = 4, align = "center",
                                 br(),
                                 
                                 #Continue with saved data
                                 actionBttn(inputId = "continue", 
                                            label = "Continue with saved data",
                                            style = "simple",
                                            color = "warning",
                                            icon = icon("fas fa-sign-in-alt"))
                          )
                          
                        ),
                        
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
                        br()
                        
                        
                        #********************************************************#
                        
               ),
               
               
               
               
               
               ###################################################################
               #  Outputs
               ###################################################################
               
               tabPanel("Output", value = "output_panel", 
                        icon = icon("fas fa-layer-group"),
                        
                        navlistPanel(id = "tabs_output",
                                     tabPanel("Table", value = "PRS_table",
                                              dataTableOutput("PRS_table")),
                                     
                                     tabPanel("Heatmap", value = "PRS_heatmap")
                        )
               )
    ) #navbarpage
  ) # fluidpage
) # taglist
