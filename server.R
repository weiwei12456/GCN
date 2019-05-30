options(shiny.maxRequestSize=30*1024^2)
library(shiny)
library(tseries)
library(vars)
library(forecast)
shinyServer(function(input, output) {
  
  output$info <- renderTable ({
    
    inFile <- input$file1
    if (is.null(inFile))
      return("No file uploaded")
    
    n <- as.numeric(input$n)
    
    B <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote,row.names = 1)
    d <- matrix(ncol = ncol(B), nrow = 1)
    colnames(d) <- colnames(B)
    for (i in 1:ncol(B)) {
      adf.test(B[,i])
      adf.res <- adf.test(B[,i])
      p <- adf.res$p.value
      v = 0
      while (p >= 0.05 & v < n) {
        dif <- diff(B[,i])
        adf.test(dif)
        adf.res <- adf.test(dif)
        p <- adf.res$p.value
        v = v+1
      } 
      if (v < n) {
        d[,i] <- v  
      } else
        d[,i] <- NA
    }
    d
  })
})
