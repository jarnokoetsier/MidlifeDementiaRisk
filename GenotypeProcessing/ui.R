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
    
    navbarPage(title = "I PRS Multi-Trait", id = "navbar",
               
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
                                                  label = "Click here to select data file",
                                                  title = "Please select a file:",
                                                  multiple = FALSE,
                                                  icon = icon("fas fa-mouse-pointer")),
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
                                            icon = icon("fas fa-mouse-pointer")),
                                 
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
                                 )
                          )
                          
                        ),
                        
                        fluidRow(
                          column(4, offset = 4, align = "center",
                                 style = "background-color:#FFFFFF;",
                                 hr(),
                                 column(4, align = "center",
                                       img(src = "MHENS_logo.png", width = "100%")),
                                 column(4, align = "center",
                                        img(src = "UM_logo.png", width = "100%")),
                                 column(4, align = "center",
                                        img(src = "EXETER_logo.png", width = "100%")),
                                 
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
                                              h1(strong("PRS Data Table")),
                                              hr(),
                                              dataTableOutput("PRS_table"),
                                              downloadButton("downloadPRS", 
                                                             "Download Table")),
                                     
                                     tabPanel("Heatmap", value = "PRS_heatmap",
                                              fluidRow(
                                                column(3,
                                                       selectInput(inputId = "cor_method",
                                                                   label = "Correlation Method",
                                                                   choices = c("pearson",
                                                                               "spearman",
                                                                               "kendall"),
                                                                   selected = "pearson")
                                                       ),
                                                column(3,
                                                       selectInput(inputId = "link_method",
                                                                   label = "Linkage Method",
                                                                   choices = c("ward.D2",
                                                                               "ward.D",
                                                                               "average",
                                                                               "complete",
                                                                               "centroid",
                                                                               "median", 
                                                                               "mcquitty"),
                                                                   selected = "ward.D2")
                                                ),
                                                
                                              ),
                                              hr(),
                                              plotOutput("heatmap_plot",
                                                         width = 1000,
                                                         height = 800,
                                                         click = NULL)%>% 
                                                withSpinner(color="red")
                                              ),
                                     tabPanel("Correlations", value = "PRS_correlation",
                                              fluidRow(
                                                column(3,
                                                       selectInput(inputId = "cor_method2",
                                                                   label = "Correlation Method",
                                                                   choices = c("pearson",
                                                                               "spearman",
                                                                               "kendall"),
                                                                   selected = "pearson")
                                                ),
                                                column(3,
                                                       uiOutput("ui_X")
                                                ),
                                                column(3,
                                                       uiOutput("ui_Y")
                                                )
                                                
                                              ),
                                              hr(),
                                              plotOutput("correlation_plot",
                                                         width = 1000,
                                                         height = 600,
                                                         click = NULL)%>% 
                                                withSpinner(color="red")
                                     )
                                     
                                     
                                     
                                     
                                     
                        )
               )
    ) #navbarpage
  ) # fluidpage
) # taglist
