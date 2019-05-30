#ui.R
library(shiny)
shinyUI(pageWithSidebar(
  headerPanel("CSV Viewer"),
  sidebarPanel(
    fileInput('file1', 'Choose CSV File',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    tags$hr(),
    checkboxInput('header', 'Header', TRUE),
    radioButtons('sep', 'Separator',
                 c(Comma=',',
                   Semicolon=';',
                   Tab='\t'),
                 'Comma'),
    radioButtons('quote', 'Quote',
                 c(None='',
                   'Double Quote'='"',
                   'Single Quote'="'"),
                 'Double Quote'),
    sliderInput("ordervalue", 
                "Max order:", 
                value = 5,
                min = 1, 
                max = 10, step = 1),
    sliderInput("at", 
                "abundance threshold", 
                value = 0.05,
                min = 0, 
                max = 1, step = 0.01),
    selectInput("dataset", "Choose a dataset:", 
                choices = c("matrixP","network","Bonferroninetwork","networkanalysis","Bonferroninetworkanalysis")),
    downloadButton('downloadData', 'Download')
  ),
  mainPanel(
  tabsetPanel(
    
    tabPanel("GCN", tableOutput('content1')),
    tabPanel("Bonferroninetwork", tableOutput('content2')),
    tabPanel("networkanalysis", tableOutput('content3')),
    tabPanel("Bonferroninetworkanalysis", tableOutput('content4'))
  )
)
))


