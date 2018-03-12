library(ggplot2)
library(plotly)
library(shiny)
library(DT)
library(edgeR)
library(statmod)
library(GO.db)
library(org.Mm.eg.db)
library(ggthemr)

ggthemr('fresh')

x <- readRDS("./Data/AAseqDEobj.RData")
group <- factor(c(rep(1:2, 4), rep(3:4, 4), rep(5:6, 4)),
                labels = c("KO_0h", "WT_0h", "KO_3h", "WT_3h", "KO_6h", "WT_6h"))
designMat <- model.matrix(~0+group)
colnames(designMat) <- levels(group)
fit <- glmQLFit(x, designMat, robust = T)
rpkmMat <- readRDS("Data/AAseqRPKMmat.RData")
rpkmMat$Time <- rep(c(0, 3, 6), each = 8)
rpkmMat$Genotype <- rep(c("KO", "WT"), 12)

function(input, output, session) {
  
  makeDEtab <- reactive({
    
    comparison <- makeContrasts(contrasts = paste(as.character(input$compChoose1),
                                                  as.character(input$compChoose2),
                                                  sep = " - "), levels = designMat)
    tr <- glmTreat(fit, contrast = comparison, 
                   lfc = log2(1.2), null = "worst.case")
    trTags <- topTags(tr, n = nrow(tr$table))
    trTab <- trTags$table
    trTab$Significant <- F
    trTab$Significant[trTab$FDR<0.05 & abs(trTab$logFC) > 1.2 ] <- T
    return(trTab)
    
  })
  
  
  makePathway <- reactive({
    
    if (input$returnPath == T){
    
    comparison <- makeContrasts(contrasts = paste(as.character(input$compChoose1),
                                                  as.character(input$compChoose2),
                                                  sep = " - "), levels = designMat)
    tr <- glmTreat(fit, contrast = comparison, 
                   lfc = log2(1.2), null = "worst.case")
    goann <- goana(tr, species = "Mm")
    keggan <- kegga(tr, species = "Mm")
    
    if (input$pathCheck == T){
      pathOut <- topKEGG(keggan, n = 30)
    } else {
      pathOut <- topGO(goann, n = 30)
    }
    
    return(pathOut)
    }
    
  })
  
  
  output$volcanoOut <- renderPlotly({
    
    if(input$scaleCheck == T){
    
      p <- ggplot(data = makeDEtab(), 
                  aes(x = logFC, y = -log10(PValue), 
                      color = Significant,
                      text = paste("Gene:", make.names(Symbol)),
                      key = make.names(Symbol))) +
        geom_point(alpha = (1/4), aes(size = logCPM)) +
        labs(x = "Log2 Fold Change", y = "-log10 q Value") +
        guides(color = guide_legend(title = "Significant")) +
        scale_color_manual(values = c("blue", "red"))
      
      p <- ggplotly(p, source = "Volcano")
      
      print(p)
        
    }else{
    
    p <- ggplot(data = makeDEtab(), 
                aes(x = logFC, y = -log10(PValue), 
                    color = Significant,
                    text = paste("Gene:", make.names(Symbol)),
                    key = make.names(Symbol))) +
      geom_point(alpha = (1/2)) +
      labs(x = "Log2 Fold Change", y = "-log10 q Value") +
      guides(color = guide_legend(title = "Significant")) +
      scale_color_manual(values = c("blue", "red"))
    
    p <- ggplotly(p, source = "Volcano")
    
    print(p)
    }
    
  })
  
  output$diffTable <- DT::renderDataTable(DT::datatable({
    
    makeDEtab()
    
  }))
  
  output$barPlotNaive <- renderPlotly({
    
    event.data <- event_data("plotly_click", source = "Volcano")
    if (is.null(event.data)) {
      
      p <- ggplot(rpkmMat, aes(x = group, y = Cxcl1)) + 
        geom_boxplot()
      p <- ggplotly(p)
      
    } else {
      
      p <- ggplot(rpkmMat, aes_string(x = "group", y = event.data$key)) + 
        geom_boxplot()
      p <- ggplotly(p)
      
    }
    
  })
  
  output$barPlotNaive2 <- renderPlotly({
    
    p <- ggplot(rpkmMat, aes_string(x = "group", y = input$geneChoose)) + 
      geom_boxplot()
    p <- ggplotly(p)
    
  })
  
  output$pathPrint <- renderPrint({
    
    makePathway()
    
  })
  
  output$TCplot <- renderPlot({
    
    q <- ggplot(data = rpkmMat, aes_string(x = "Time", y = input$TCseq, color = "Genotype")) +
      geom_point() +
      stat_smooth(method = 'lm', formula = y ~ poly(x, 2), size = 1, alpha = 0.1)
    
    print(q)
    
  })

  output$dataDownload <- DT::renderDataTable(t(rpkmMat))
  output$dataSelect <- renderPrint(input$dataDownload_rows_selected)
    
  RPKMdown <- t(rpkmMat)
  selectedRow <- eventReactive(input$dataDownload_rows_selected,{
    row.names(RPKMdown)[c(input$dataDownload_rows_selected)]
  })
  
  output$selected <- renderText({
    selectedRow()
  })
  
  output$downloadData <- downloadHandler(
    
    filename = "seqDataDownload.csv",
    content = function(file) {
      write.csv(RPKMdown[row.names(RPKMdown)[c(input$dataDownload_rows_selected)],],
                file)
#      write.csv(data.frame(RPKMdown[row.names(RPKMdown) %in% row.names(RPKMdown)[c(input$dataDownload_rows_selected)]]),
#                  file)
      
    })
  
#  getIndex <- observeEvent(input$getData,{
 #   indexVec <- c(input$dataDownload_rows_selected)
  #  RPKMdown <- RPKMdown[indexVec,]
  #)
  
  #output$checkDown <- renderPrint(length(RPKMdown))
  

  }