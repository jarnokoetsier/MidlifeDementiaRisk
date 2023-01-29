Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

server <- function(input, output, session){
  options(shiny.maxRequestSize=1000000*1024^2)
  
  ##############################################################################
  #Data selection
  ##############################################################################
  
  
  #Hide tabs
  hideTab("navbar", target = "output_panel")
  
  #Info panel
  observeEvent(input$infopanel1, {
    sendSweetAlert(
      session = session,
      title = "Information",
      text = "Information should be provided here",
      type = "info"
    )
  })
  
  
  #****************************************************************************#
  #Saved data
  #****************************************************************************#
  
  #Enter ID to retrieve saved data
  observeEvent(input$continue, {
    inputSweetAlert(
      session = session,
      inputId = "idcontinue",
      title = "Continue",
      text = "Enter the cohort name (i.e., the name of the PLINK file without .bim/.bed/.fam extension) 
      to continue.",
      input = "text",
      inputPlaceholder = "Cohort name",
      btn_labels = c("Continue", "Cancel")
    )
    
  })
  
  volumes = getVolumes()
  observe({
    shinyFileChoose(input, "bfile", roots = volumes, session = session)
  })
  
  #bfile
  bfile <- eventReactive(input$startAnalysis,{
    req(input$bfile)
    bfile <- as.character(parseFilePaths(volumes, input$bfile)$datapath)
    return(bfile)
  })
  
  continue_text <- reactive({
    req(input$idcontinue)
    return(input$idcontinue)
  })
  
  # Select traits
  traits_selected <- eventReactive(input$startAnalysis,{
    if (input$all_traits == TRUE){
      traits_selected <- Traits
    }
    if (input$all_traits == FALSE){
      traits_selected <- input$traits_input
    }
    
    return(traits_selected)
  })
  
  # Calculate PRS
  observeEvent(input$startAnalysis,{
   
    
    path <- reactive({
      req(traits_selected())
      req(bfile())
      path <- str_remove(bfile(), ".bed")
      path <- str_remove(path, ".bim")
      path <- str_remove(path, ".fam")
      return(path)
    })
    observe({
      req(bfile())
      print(bfile())
    })
   
    for (t in traits_selected()){
      predPRS(bfile = wslPath(path()), 
              Trait = t, 
              OverlapSNPsOnly=FALSE, 
              Force = FALSE)
    }
  })

  # Get cohort name

  cohortName <- reactive(NULL)
 observeEvent(input$startAnalysis,{
   cohortName <- reactive({
     if (length(bfile()) > 0){
       path <- str_remove(bfile(), ".bed")
       path <- str_remove(path, ".bim")
       path <- str_remove(path, ".fam")
       cohortName <- str_remove(path, ".*/")
       
       return(cohortName)
     }
   })
 })   

 observeEvent(input$idcontinue,{
   cohortName <- reactive({
     if (length(continue_text()) > 0){
       cohortName <- continue_text()
       return(cohortName)
     }
     
   })
 })

  

  
  
  # Get PRS results
  PRS_result <- reactive({
    req(cohortName())
    PRS_result <-  collect_all_PRS(cohort = cohortName())
    return(PRS_result)
  })
  
  # Output table
  output$PRS_table <- DT::renderDataTable({
    return(PRS_result())
  }, server=TRUE,
  options = list(pageLength = 5), rownames= FALSE)
  
  # Go the next step
  observeEvent(if (length(PRS_result()) > 0){input$idcontinue}, {
    
    # Success message
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully selected!",
      type = "success")
    
    # Go to next tab
    updateTabsetPanel(session, "navbar",
                      selected = "output_panel")
    
    # Show next tab
    showTab("navbar", target = "output_panel")
    
    
  })
  
  # Go the next step
  observeEvent(if (length(PRS_result()) > 0){input$startAnalysis}, {
    
    # Success message
    sendSweetAlert(
      session = session,
      title = "Success!",
      text = "Data successfully selected!",
      type = "success")
    
    # Go to next tab
    updateTabsetPanel(session, "navbar",
                      selected = "output_panel")
    
    # Show next tab
    showTab("navbar", target = "output_panel")
    
  })
}