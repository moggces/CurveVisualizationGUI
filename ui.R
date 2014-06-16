
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("General Concentration-Response Data Plotting"),
  
  # Sidebar with a slider input for number of observations
  sidebarPanel(
    h4('Mode'),
    radioButtons("mode", "Select a pathway display mode:",
                 choices = list("parallel"="parallel",
                                "overlay"="overlay", 
                                "parallel pathway + overlay cmpd"="mixed")),
    tags$br(),
    
    h4('Input'),
    fileInput('file1', 'Import conc.-resp. data', multiple=TRUE),
    
    h4('Compound filter'),
    tags$textarea(id="cmpds", rows=3, cols=1, ""),
    
    tags$hr(),
    
    sliderInput("widthpx", 
                "width pixel/colmum", min = 50, max = 1000, value = 300, step=50),
    sliderInput("heightpx", 
                "height pixel/colmum", min = 50, max = 1000, value = 300, step=50),
    
    tags$br(),
    
    h4('Pathway'),
    wellPanel (
      uiOutput("pathways")
    ),
    
    h4('Pathway readout options'),
    wellPanel (
      uiOutput("options"),
      #checkboxInput("isOneAssay", "multiplex cytotoxicity as one assay", TRUE)
      checkboxInput("useParent", "use parent tag", FALSE)
    ),
    
    h4('Curve plotting options'),
    wellPanel (
      uiOutput("plot_options"),
      checkboxInput("showOutlier", "cross outliers", TRUE)
    ),
    br(),
    downloadButton('downloadPlot', 'Save Plot')
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    
    tabsetPanel(
      tabPanel("Text", dataTableOutput('temp')),
      tabPanel( 'Data', dataTableOutput('contents')),
      tabPanel( "Plot", plotOutput("plot", height="auto", width="500%"))
      
      #tabPanelAbout()
    )
  )
))
