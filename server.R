
library(shiny)
library(plyr)
library(reshape2)
library(ggplot2)

source("./source/load.R",  local=TRUE)
source("./source/get.R",  local=TRUE)
options(shiny.maxRequestSize=30*1024^2)

# to add a new analysis method, edit plot_options, get_melt_data(), get_plot(), get_blank_data()
# know bugs: if Tox21_110433_1 alone ; there will be error 
################

shinyServer(function(input, output) {
   
  data_loader <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return (NULL)
    result <- NULL
    
   
      result <- lapply(1:nrow(inFile),
                       function (x)
                       {
                         path <- inFile[[x, 'datapath']]
                         load_input_file(path)
                       }
      )
      result <- do.call("rbind.fill", result)

    
#     path <- inFile$datapath
#     #result <- data.frame()
#     result <- NULL
#     result <- load_input_file(path)
    return(result)
  })

  
  data_chemical <- reactive({
    result <- NULL
    ids <- input$cmpds
    
    #if (is.null(ids)) return (NULL)
    
    result <- unlist(strsplit(ids, "\n", perl=TRUE))
    return(result)
  })
  
  output$pathways <- renderUI({
    result <- data_loader()
    if (is.null(result)) return(NULL)
    pp <- list()
    pathway <- unique(result$pathway)
    pp  <-   split(pathway, 1:length(pathway))
    names(pp) <- pathway
    selectInput("paths", "Select pathways to show:", 
                choices  = pp, 
                multiple = TRUE)
  })
  
  output$options <- renderUI({
    result <- data_loader()
    if (is.null(result)) return(NULL)
    pp <- list()
    readout <- unique(result$readout)
    pp  <-   split(readout, 1:length(readout))
    names(pp) <- readout
    
    selectInput("opts", "Select readout options:", 
                choices  = pp,
                multiple = TRUE)
  })

  output$plot_options <- renderUI({
    selectInput("plt_opts", "Select line plotting options:", 
                choices  = list( "raw"="raw", "curvep"="curvep",
                                 "Hill 4-point"="hill",
                                 "Hill 4-point(fred)"="hill_fred"),
                selected=c("raw"),
                multiple = TRUE)
  })

  data_filter <- reactive ({
    
    result <- data_loader()
    chemicals <- data_chemical()
    
    options <- input$opts
    pathways <- input$paths
  
    readout_pat <- paste(options, sep="", collapse="|")
    pathway_pat <- paste(pathways, sep="", collapse="|")
    
    if (readout_pat == '' ) { result <- data.frame()  }  
    if (pathway_pat == '')  {result <- data.frame() } 
    #if (sum(chemicals == '') > 0 ) result <- data.frame()
    if ( is.null(chemicals) ) { result <- data.frame() } ### not working; don't know why
    
    #result <- result[grep(pathway_pat, result$pathway, perl=TRUE),]
    result <- result[grepl(pathway_pat, result$pathway, perl=TRUE) | grepl(pathway_pat, result$parent, perl=TRUE),]
    result <- result[grep(readout_pat, result$readout, perl=TRUE),]
    result <- get_qhts_data(chemicals, result)
    
    return(result)
  })
  
  data_melter <- reactive ({
    mode <- input$mode
    plot_options <- input$plt_opts
    show_outlier <- input$showOutlier
    if (show_outlier) plot_options <- c(plot_options, 'mask')
    use_parent <- input$useParent
    
    qhts <- data_filter()
    
    if (use_parent) qhts <- get_long_format(qhts)
    
    result <- get_melt_data(qhts, resp_type=unique(c('raw', plot_options)))
    #result$display_name <- paste(result$CAS, "|\n", result$Chemical.ID, sep="")
    result$display_name <- paste(result$Chemical.ID, "|\n", result$Chemical.Name, sep="")
    result <- result[order(result$display_name),]
    return(result)
  })
  
#   getVarHeight <- function() {
#     qhts <- data_filter()
#     nrow <- length(unique(qhts$Chemical.ID))
#     mode <- input$mode
#     if (mode == 'overlay')
#     {
#       return(nrow * 150)
#     } else if (mode == 'parallel')
#     {
#       return(nrow * 300)
#     }
#   }
#   
#   getVarWidth <- function() {
#     qhts <- data_filter()
#     mode <- input$mode
#     ncol <- length(unique(qhts$pathway))
#     if (mode == 'overlay')
#     {
#       #return("auto")
#       return(ncol * 300)
#     } else if (mode == 'parallel')
#     {
#       return(ncol * 600)
#     }
#   }
  
  getVarHeight <- reactive({
    qhts <- data_filter()
    nrow <- length(unique(qhts$Chemical.ID))
    mode <- input$mode
    heightpx <- input$heightpx
    if (mode == 'overlay' )
    {
      return(nrow * heightpx) # 350
    } else if (mode == 'parallel' | mode == 'mixed')
    {
      return(nrow * heightpx) # 300
    }
  })
  
  getVarWidth <- reactive({
    qhts <- data_filter()
    mode <- input$mode
    ncol <- length(unique(qhts$pathway))
    widthpx <- input$widthpx
    if (mode == 'overlay')
    {
      #return("auto")
      return(1000)
    } else if (mode == 'parallel'  | mode == 'mixed')
    {
      return(ncol * widthpx)
    }
  })
  
  output$plot <- renderPlot({
    mode <- input$mode
    plot_options <- input$plt_opts
    show_outlier <- input$showOutlier
    if (show_outlier) plot_options <- c(plot_options, 'mask')
    
    result <- data_melter()
    p <- get_plot(result, mode=mode, plot_options=plot_options, fontsize=20, pointsize=3)
    
    if (mode == 'overlay')
    {
      p <- p  + theme_bw(base_size = 20) + facet_wrap(~ display_name  , ncol=2)
    } else if (mode == 'parallel')
    {
      p <- p + theme_bw(base_size = 20) + facet_grid(display_name ~ pathway)
    } else if (mode == 'mixed')
    {
      p <- p  + theme_bw(base_size = 20) + facet_wrap(~ pathway  , ncol=2)
    }
    print(p)
    #print(select_plot())
  }, height=getVarHeight, width=getVarWidth)

  output$downloadPlot <- downloadHandler(
    filename = function() { paste(as.numeric(as.POSIXct(Sys.time())), ".pdf", sep="") },
    content = function(file) {
      result <- data_melter()
      mode <- input$mode
      plot_options <- input$plt_opts
      show_outlier <- input$showOutlier
      if (show_outlier) plot_options <- c(plot_options, 'mask')
      
      #height <- getVarHeight()
      #width <- getVarWidth()
      #pdf(file=file, width=width/96, height=height/96)
      
      #pdf(file=file, paper = "USr")
      
      pdf(file=file)
      
      if (mode == 'overlay')
      {
        n_page <- 6
        nn <- unique(result$display_name)
        pages <- split(nn, ceiling(seq_along(nn)/n_page))
        lapply(names(pages), function (x) {
          sub <- result[result$display_name %in% pages[[x]],]
          p <- get_plot(sub, mode=mode, plot_options=plot_options, fontsize=8, pointsize=1)
          p <- p  + theme_bw(base_size = 8) + facet_wrap(~ display_name  , ncol=2, nrow=3) 
          print(p)
        })
        
      } else if (mode == 'parallel')
      {
        n_page <- 6
        result <- get_blank_data(result, n_page)
        nn <- unique(result$display_name)
        pages <- split(nn, ceiling(seq_along(nn)/n_page))
        lapply(names(pages), function (x) {
          sub <- result[result$display_name %in% pages[[x]],]
          p <- get_plot(sub, mode=mode, plot_options=plot_options, fontsize=8, pointsize=1)
          p <- p + theme_bw(base_size = 8) + facet_grid(display_name ~ pathway)
          print(p)
        })
        
      } else if (mode == 'mixed') ## not done
      {
        
      }
      
      dev.off()
    }
  )

  output$contents <- renderDataTable({
    result <- data_filter()
    #result <- data_loader()
    return(result)
  })
  
  output$temp <- renderDataTable({
    result <- data_loader()
    #result <- data_filter()
    #result <- data_melter()
    #result <- unique(result$parent)
    return(result)
  })
  
})
