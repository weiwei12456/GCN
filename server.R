options(shiny.maxRequestSize=30*1024^2)
shinyServer(function(input, output) {
  d1 <- reactive({
    
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    dt0 <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, row.names = 1)
    parameter <- data.frame(
      Name = c("ordervalue", 
               "at"),
      Value = as.character(c(input$ordervalue, 
                             input$at)), 
      stringsAsFactors=FALSE)
    pam <- as.matrix(parameter)
    ordervalue <- as.numeric(input$ordervalue)
    at <- as.numeric(input$at)
    B <- dt0 # check
    head(B)
    # confirm d of stationary
    library(tseries)
    library(vars)
    library(forecast)
    d <- matrix(ncol = ncol(B),nrow = 1)
    colnames(d) <- colnames(B)
    for (i in 1:ncol(B)) {
      adf.test(B[,i])
      adf.res <- adf.test(B[,i])
      p <- adf.res$p.value
      v = 0
      while (p >= 0.05 & v < 5) {
        dif <- diff(B[,i])
        adf.test(dif)
        adf.res <- adf.test(dif)
        p <- adf.res$p.value
        v = v+1
      } 
      if (v < 5) {
        d[,i] <- v  
      } else
        d[,i] <- 0
    }
    d
    #stationary 
    a <- matrix(0,nrow = (nrow(B)-max(d)),ncol = ncol(B))#Build the null matrix for accepting value from the checking results.
    ncol(B)
    for  (i in 1:ncol(B)){
      dif <- B[,i]
      v = 0
      n = d[1,i] #####################
      adf.test(dif)
      adf.res <- adf.test(dif)
      p <- adf.res$p.value
      while (v < n){
        dif <- diff(dif)
        v <- v+1
        adf.test(dif)
        adf.res <- adf.test(dif)
        p <- adf.res$p.value
      } 
      if (p < 0.05) {
        mdif<- as.matrix(dif)
        a[,i] <- mdif[1:(nrow(B)-max(d)),]
      }
    } 
    a# stationary matrix
    
    #stationary check
    colnames(a) <- colnames(B)# name the original rowname and colname
    checkunstaionary <- a[,which(colSums(a == 0)==NROW(a))]#take out unstationary sequences
    checkunstaionary
    dim(checkunstaionary)
    a <- a[,which(colSums(a == 0) < NROW(a))]#take out unstationary sequences
    
    #Define the lag order according to AIC
    library(vars)
    lagvector <- matrix(0, nrow = ncol(a), ncol = 1)
    rownames(lagvector) <- colnames(a)
    for (i in 1:ncol(a)) {
      var<-VARselect(a[,i],lag.max = ordervalue)
      vars<- as.numeric(var$selection) 
      lagvector[i] <- vars[1]
    }
    
    
    b = t(a)
    
    mat1 <- matrix(3,nrow = nrow(b),ncol = ncol(a))
    matp1 <- matrix(NA,nrow = nrow(b),ncol = ncol(a))
    
    CasualnetworkP1 <- matrix(0, nrow = nrow(b)*nrow(b), ncol = 4)
    
    colnames(CasualnetworkP1) <- letters[1:4] 
    colnames(CasualnetworkP1) <- c("OTUT","P","OTUS","Inter")
    
    
    # create null matrix to accept the values of grangertest.
    for(i in 1:ncol(a)){
      for(j in 1:nrow(b)){
        if(i ==j  ){
          mat1[i,j] <- NA
        }else{
          x = a[,i]
          y = b[j,]
          stat_data <- data.frame(x, y)
          ts.data <-ts(stat_data,frequency=1,start=c(1))
          ts.data
          library(lmtest)
          fit0 <- try({grangervalue <- grangertest(x ~ y,order = lagvector[j,], data = ts.data)
          P0 <- grangervalue$P}, silent = TRUE)
          if ('try-error' %in% class(fit0)) {
            P0 <- matrix(5,nrow = 2, ncol = 1) 
          }else {}
          
          P1 <- matrix(P0)[2,1]
          matp1[i,j] <- P1
          
          
          if (P1 > at ){  
            mat1[i,j] <- NA 
          }else  
            mat1[i,j] <- P1 
        }
        
      }
    }
    mat1
    matp1
    
    # This matrix was already an adjacency matrix.
    rownames(mat1) <- rownames(b)
    colnames(mat1) <- colnames(a)
    rownames(matp1) <- rownames(b)
    colnames(matp1) <- colnames(a)
    mat1
  })
  
  d2 <- reactive ({
     at <- as.numeric(input$at)
     
    trans <- d1()
    CasualnetworkP1 <- matrix(0, nrow = nrow(trans)*nrow(trans), ncol = 4)
    trans[is.na(trans)]<-100
    
    colnames(CasualnetworkP1) <- c("OTUT","P","OTUS","Inter")
    
    for(j in 0:(nrow(trans)-1)){
      namematrix <- matrix(rownames(trans)[j+1],nrow = nrow(trans),ncol = ncol(trans))
      CasualnetworkP1[(1+(j*nrow(trans))):((j+1)*nrow(trans)),1]<-namematrix[,1]
      CasualnetworkP1[(1+(j*nrow(trans))):((j+1)*nrow(trans)),3]<-colnames(trans) 
      CasualnetworkP1[(1+(j*nrow(trans))):((j+1)*nrow(trans)),2]<-trans[j+1,]
    }
    
    for (i in 1:nrow(CasualnetworkP1)) {
      if (as.numeric(CasualnetworkP1[i,2]) < at) {
        CasualnetworkP1[i,4] <- as.numeric(CasualnetworkP1[i,2])
      } else 
        CasualnetworkP1[i,4] <- NA
    }
    
    checkmatrix<- matrix(nrow = nrow(trans),ncol = ncol(trans))
    
    for (e in nrow(trans)) {
      for (f in ncol(trans)) {
        if (trans[e,f]==5)  {
          checkmatrix[e,f] <- trans[e,f]
        }
      }
    }#check value of mat1, 5 means the calculation halted from here
    checkmatrix
    networkfileat <- CasualnetworkP1[which(as.numeric(CasualnetworkP1[,4]) != 0),]
  })
  d3 <- reactive ({
    at <- as.numeric(input$at)
    trans <- d1()
    networkfileat <- d2()
    Bocutoff <- at/choose(ncol(trans),2)
    networkfilebo <- networkfileat[which(as.numeric(networkfileat[,4]) < Bocutoff),]
    networkfilebo
  })  
  d4 <- reactive ({
    
    networkfilebolag5 <- d2()
    a <- d1()
    otutotal <- matrix(0, ncol = 9, nrow = ncol(a))
    otutotal[,1] = colnames(a)
    networkfilebolag5
    for (i in 1:nrow(otutotal)) {
      d<-otutotal[i,1]
      otuse1 = as.matrix(networkfilebolag5[which(networkfilebolag5[,1] %in% d),],ncol = 4)
      otuse2 = as.matrix(networkfilebolag5[which(networkfilebolag5[,3] %in% d),],ncol = 4)
      replicate = matrix(otuse1[which(otuse1[,3] %in% otuse2[,1]),],ncol = 4)
      dim(replicate)
      difference <- dim(otuse1)+dim(otuse2)-dim(replicate)
      otutotal[i,4] <- difference[1]
      otutotal[i,2] <- dim(otuse1)[1]
      otutotal[i,3] <- dim(otuse2)[1]
      
    }
    colnames(otutotal) <- c("OTU","indegree","outdegree","neighbors","cs","crecip","no","ni", "ncf")
    otutotal[,5] <- as.numeric(otutotal[,2])/(as.numeric(otutotal[,3])+as.numeric(otutotal[,2]))
    otutotal[,6] <- (as.numeric(otutotal[,3])+as.numeric(otutotal[,2])-as.numeric(otutotal[,4]))/(as.numeric(otutotal[,3])+as.numeric(otutotal[,2]))
    
    otutotal[,7] <- (as.numeric(otutotal[,4])-as.numeric(otutotal[,2]))
    otutotal[,8] <- (as.numeric(otutotal[,4])-as.numeric(otutotal[,3]))
    otutotal[,9] <- as.numeric(otutotal[,7])-as.numeric(otutotal[,8])
    otutotal 
  })  
  d5 <- reactive ({
    
    networkfilebolag5 <- d3()
    a <- d1()
    otutotal <- matrix(0, ncol = 9, nrow = ncol(a))
    otutotal[,1] = colnames(a)
    networkfilebolag5
    
    for (i in 1:nrow(otutotal)) {
      d<-otutotal[i,1]
      otuse1 = as.matrix(networkfilebolag5[which(networkfilebolag5[,1] %in% d),],ncol = 4)
      otuse2 = as.matrix(networkfilebolag5[which(networkfilebolag5[,3] %in% d),],ncol = 4)
      replicate = matrix(otuse1[which(otuse1[,3] %in% otuse2[,1]),],ncol = 4)
      dim(replicate)
      difference <- dim(otuse1)+dim(otuse2)-dim(replicate)
      otutotal[i,4] <- difference[1]
      otutotal[i,2] <- dim(otuse1)[1]
      otutotal[i,3] <- dim(otuse2)[1]
      
    }
    colnames(otutotal) <- c("OTU","indegree","outdegree","neighbors","cs","crecip","no","ni", "ncf")
    otutotal[,5] <- as.numeric(otutotal[,2])/(as.numeric(otutotal[,3])+as.numeric(otutotal[,2]))
    otutotal[,6] <- (as.numeric(otutotal[,3])+as.numeric(otutotal[,2])-as.numeric(otutotal[,4]))/(as.numeric(otutotal[,3])+as.numeric(otutotal[,2]))
    
    otutotal[,7] <- (as.numeric(otutotal[,4])-as.numeric(otutotal[,2]))
    otutotal[,8] <- (as.numeric(otutotal[,4])-as.numeric(otutotal[,3]))
    otutotal[,9] <- as.numeric(otutotal[,7])-as.numeric(otutotal[,8])
    otutotal 
  }) 
 
  output$content1 <- renderTable({
    if(is.null(input$file1))
    {
      return("No file uploaded");
    }else
    {
      d2()
    }
  })
  output$content2 <- renderTable({
    if(is.null(input$file1))
    {
      return("No file uploaded");
    }else
    {
      d3()
    }
  })
  output$content3 <- renderTable({
    if(is.null(input$file1))
    {
      return("No file uploaded");
    }else
    {
      d4()
    }
  })
  output$content4 <- renderTable({
    if(is.null(input$file1))
    {
      return("No file uploaded");
    }else
    {
      d5()
    }
  })
  datasetInput <- reactive({
    switch(input$dataset,
           "network" = d2(),
           "matrixP" = d1(),
           "Bonferroninetwork" = d3(),
           "networkanalysis" = d4(),
           "Bonferroninetworkanalysis" = d5()
    )
  })
  output$downloadData <- downloadHandler(
    
    filename = function() { paste(input$dataset, '.csv', sep='') },
    content = function(file) {
      
      write.csv(datasetInput(), file)
    }
  )
  
})
