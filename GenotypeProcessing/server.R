Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

server <- function(input, output, session){
  options(shiny.maxRequestSize=1000000*1024^2)
  
  #****************************************************************************#
  # Data selection
  #****************************************************************************#
  
  
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
  
  volumes = getVolumes()
  observe({
    shinyFileChoose(input, "bfile", roots = volumes, session = session)
  })
  
  
  #****************************************************************************#
  # Continue with saved data
  #****************************************************************************#
  
  
  observeEvent(input$continue, {
    
    # Enter cohort name
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
    
    # Cohort name
    cohortName <- reactive({
      req(input$idcontinue)
      return(input$idcontinue)
    })
    
    # Get PRS results
    PRS_result <- reactive({
      req(cohortName())
      PRS_result <-  collect_all_PRS(cohort = cohortName())
      return(PRS_result)
    })
    
    
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
    
    # Output table
    output$PRS_table <- DT::renderDataTable({
      return(PRS_result())
    }, server=TRUE,
    options = list(pageLength = 10), rownames= TRUE)
    
    
    #download table of meta data
    output$downloadPRS <- downloadHandler(
      filename = "PRS_results.csv",
      content = function(file){
        write.csv(PRS_result(), file)
      }
    )
    
    # Heatmap
    output$heatmap_plot <- renderPlot({
      req(input$cor_method)
      req(input$link_method)
      
      # Calculate correlations
      corrDF <- as.data.frame(correlate(PRS_result(), diagonal = 1, method = input$cor_method))
      rownames(corrDF) <- corrDF$term
      corrDF <- corrDF[,-1]
      
      # Format data for plotting
      plotCor <- gather(corrDF)
      plotCor$key1 <- rep(rownames(corrDF), ncol(corrDF))
      
      # Perform clustering to get sample order
      model <- hclust(as.dist(1-abs(corrDF)), input$link_method)
      order <- model$labels[model$order]
      plotCor$key <- factor(plotCor$key, levels = order)
      plotCor$key1 <- factor(plotCor$key1, levels = order)
      
      # Main correlation plot
      main <- ggplot(plotCor) +
        geom_tile(aes(x = key, y = key1, fill = value)) +
        scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                             limits = c(-1,1)) +
        xlab(NULL) +
        ylab(NULL) +
        labs(fill = paste0(str_to_title(input$cor_method),"\nCorrelation"))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15),
              axis.text.y = element_text(size = 15))
      
      # Dendrogram
      dendroPlot <- ggplot() +
        geom_tile(data = plotCor, aes(x = as.numeric(key), y = 1), fill = "white", alpha = 0) +
        geom_dendro(model,xlim = c(1,length(unique(plotCor$key)))) +
        theme_void() +
        theme(legend.position = "none") 
      
      # Combine plots
      p <- dendroPlot + main +
        plot_layout(nrow = 2, ncol = 1,
                    heights = c(1,6))
      
      # Return output
      return(p)
      
    })
    
    #x-axis
    output$ui_X <- renderUI({
      selectInput(inputId = "X",
                  label = "X-axis",
                  choices = colnames(PRS_result()),
                  selected = colnames(PRS_result())[1])
    })
    
    # y-axis
    output$ui_Y <- renderUI({
        selectInput(inputId = "Y",
                    label = "Y-axis",
                    choices = colnames(PRS_result()),
                    selected = colnames(PRS_result())[2])
      
    })
    
    # Correlation plot
    output$correlation_plot <- renderPlot({
      req(input$cor_method2)
      req(input$X)
      req(input$Y)
      
      plotScatter <- data.frame(X = PRS_result()[,input$X],
                                Y = PRS_result()[,input$Y])
      
      corrValue = cor(plotScatter$X,plotScatter$Y, method = input$cor_method2,
                      use = "pairwise.complete.obs")
      pValue = cor.test(plotScatter$X,plotScatter$Y, method = input$cor_method2,
                        use = "pairwise.complete.obs")$p.value
      p <- ggplot(plotScatter) +
        geom_point(aes(x = X, y = Y), color = "#DC3535") +
        xlab(input$X) +
        ylab(input$Y) +
        ggtitle(paste0(input$X, " vs ", input$Y),
                subtitle = paste0("Corr. coeff = ", round(corrValue,3), ", p-value = ", round(pValue,3))) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5,
                                        face = "bold",
                                        size = 24),
              plot.subtitle = element_text(hjust = 0.5,
                                           size = 18,
                                           face = "italic"),
              axis.text = element_text(size = 15),
              axis.title = element_text(size = 18))
      
      return(p)
      
      
    })
  }) #observeEvent
  

  
  #****************************************************************************#
  # Start from PLINK files
  #****************************************************************************#
  
  observeEvent(input$startAnalysis, {
    
    showModal(modalDialog(title = h4(strong("Calculating PRS..."),
                                     align = "center"), 
                          footer = NULL,
                          h5("This might take a while. Please be patient.", 
                             align = "center")))
    #bfile
    bfile <- reactive({
      req(input$bfile)
      bfile <- as.character(parseFilePaths(volumes, input$bfile)$datapath)
      return(bfile)
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
    
    # Path to data
    path <- reactive({
      req(traits_selected())
      req(bfile())
      path <- str_remove(bfile(), ".bed")
      path <- str_remove(path, ".bim")
      path <- str_remove(path, ".fam")
      return(path)
    })
    
    # Cohort name
    cohortName <- reactive({
      req(path())
      cohortName <- str_remove(path(), ".*/")
      
      return(cohortName)
    })
    
    # Get PRS results
    PRS_result <- reactive({
      req(cohortName())
      req(traits_selected())
      req(path())
      
      for (t in traits_selected()){
        predPRS(bfile = wslPath(path()), 
                Trait = t, 
                OverlapSNPsOnly=FALSE, 
                Force = FALSE)
      }
      
      PRS_result <-  collect_all_PRS(cohort = cohortName())
      removeModal()
      return(PRS_result)
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
    
    # Output table
    output$PRS_table <- DT::renderDataTable({
      return(PRS_result())
    }, server=TRUE,
    options = list(pageLength = 10), rownames= TRUE)
    
    
    #download table of meta data
    output$downloadPRS <- downloadHandler(
      filename = "PRS_results.csv",
      content = function(file){
        write.csv(PRS_result(), file)
      }
    )

    # Heatmap
    output$heatmap_plot <- renderPlot({
      req(input$cor_method)
      req(input$link_method)
      
      # Calculate correlations
      corrDF <- as.data.frame(correlate(PRS_result(), diagonal = 1, method = input$cor_method))
      rownames(corrDF) <- corrDF$term
      corrDF <- corrDF[,-1]
      
      # Format data for plotting
      plotCor <- gather(corrDF)
      plotCor$key1 <- rep(rownames(corrDF), ncol(corrDF))
      
      # Perform clustering to get sample order
      model <- hclust(as.dist(1-abs(corrDF)), input$link_method)
      order <- model$labels[model$order]
      plotCor$key <- factor(plotCor$key, levels = order)
      plotCor$key1 <- factor(plotCor$key1, levels = order)
      
      # Main correlation plot
      main <- ggplot(plotCor) +
        geom_tile(aes(x = key, y = key1, fill = value)) +
        scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                             limits = c(-1,1)) +
        xlab(NULL) +
        ylab(NULL) +
        labs(fill = paste0(str_to_title(input$cor_method),"\nCorrelation"))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 15),
              axis.text.y = element_text(size = 15))
      
      # Dendrogram
      dendroPlot <- ggplot() +
        geom_tile(data = plotCor, aes(x = as.numeric(key), y = 1), fill = "white", alpha = 0) +
        geom_dendro(model,xlim = c(1,length(unique(plotCor$key)))) +
        theme_void() +
        theme(legend.position = "none") 
      
      # Combine plots
      p <- dendroPlot + main +
        plot_layout(nrow = 2, ncol = 1,
                    heights = c(1,6))
      
      # Return output
      return(p)
      
    })
    
    #x-axis
    output$ui_X <- renderUI({
      selectInput(inputId = "X",
                  label = "X-axis",
                  choices = colnames(PRS_result()),
                  selected = colnames(PRS_result())[1])
    })
    
    # y-axis
    output$ui_Y <- renderUI({
      selectInput(inputId = "Y",
                  label = "Y-axis",
                  choices = colnames(PRS_result()),
                  selected = colnames(PRS_result())[2])
      
    })
    
    # Correlation plot
    output$correlation_plot <- renderPlot({
      req(input$cor_method2)
      req(input$X)
      req(input$Y)
      
      plotScatter <- data.frame(X = PRS_result()[,input$X],
                                Y = PRS_result()[,input$Y])
      
      corrValue = cor(plotScatter$X,plotScatter$Y, method = input$cor_method2,
                      use = "pairwise.complete.obs")
      pValue = cor.test(plotScatter$X,plotScatter$Y, method = input$cor_method2,
                        use = "pairwise.complete.obs")$p.value
      p <- ggplot(plotScatter) +
        geom_point(aes(x = X, y = Y), color = "#DC3535") +
        xlab(input$X) +
        ylab(input$Y) +
        ggtitle(paste0(input$X, " vs ", input$Y),
                subtitle = paste0("Corr. coeff = ", round(corrValue,3), ", p-value = ", round(pValue,3))) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5,
                                        face = "bold",
                                        size = 24),
              plot.subtitle = element_text(hjust = 0.5,
                                           size = 18,
                                           face = "italic"),
              axis.text = element_text(size = 15),
              axis.title = element_text(size = 18))
      
      return(p)
      
      
    })
    
  }) #observeEvent
  
}

