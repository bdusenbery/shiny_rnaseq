# server.R

# required libraries. 
library("shiny")
library("ggplot2")
library("dplyr")
library("ggvis")
library("shinyIncubator")
# data
dat <- readRDS("data/iMN_diff.rds")
fpkm <- readRDS("data/genesFpkm.rds")
fpkm$sample_good_name <- factor(fpkm$sample_good_name, levels =c(
  "Embryonic MNs", "iMNs", "ESC MNs", "iPSC MNs", "MEF1", "MEF2", "MEF3", "MEF4","ESC", "iPSC"), 
  ordered=T)


source("helper.R")
#allGeneIds <- as.character(fpkm$gene)

shinyServer(
  function(input, output, session) {
    
    # selecting the comparisons from the big data table is probably the slowest operation and should only be done once
    # until the selection is changed awesome. this helps a lot. 
    dataSelect <- reactive({
      withProgress(session, {
        setProgress(message = "Calculating, please wait", 
                    detail = "This may take a few moments...")
        setProgress(detail = "Almost there... ")
      })
      
      sampleXname <- input$xAxis
      
      sampleYname <- input$yAxis    
      
      # select samples
      datRes <- filter(dat, ((sample_1 == sampleYname & sample_2 == sampleXname) |
                               (sample_1 == sampleXname & sample_2 == sampleYname)))
            
    })
    
    dataInput <- reactive({
      
      sampleXname <- input$xAxis
       
      sampleYname <- input$yAxis      
      
      # select fold change between desired quantities and use pseudo (+1 before log2 transformation)
      if(input$pseudoFC){
        datRes <- filter(dataSelect(), (abs(pseudoFC) >= input$minFoldChange[1] ))
      }else{
        datRes <- filter(dataSelect(), (abs(log2.fold_change.) >= input$minFoldChange[1] ))
      }
      
      # if show only sig genes select desired FDR.
      if(input$sigOnly){
        datRes <- filter(datRes, (q_value <= input$minSig))
      }
      
      # select fpkm values to show, add columns with Raw FPKM values that we can print eventually, 
      # calculate proper fold change ( it is fine to use the pre-calculated
      # pseudoFC above for the inital selection as this takes the abs.  
      
      
      #
      if((datRes[1, "sample_1"] == sampleXname) & (datRes[1, "sample_2"] == sampleYname)){
        datRes <- filter(datRes, (value_1 >= input$minXfpkm &
                                    value_1 <= input$maxXfpkm &
                                    value_2 >= input$minYfpkm &
                                    value_2 <= input$maxYfpkm))%.%
          mutate(axisXprinted = value_1, 
                 axisYprinted = value_2, 
                 pseudoFC2 = log2((1+value_2)/(1+value_1)))
      }else{
        datRes <- filter(datRes, (value_1 >= input$minYfpkm &
                                    value_1 <= input$maxYfpkm & 
                                    value_2 >= input$minYfpkm &
                                    value_2 <= input$maxYfpkm))%.%
          mutate(axisYprinted = value_1,
                 axisXprinted = value_2, 
                 pseudoFC2 = log2((1+value_1)/(1+value_2)))
      }
      
      # transform fpkm values (for plotting only)
      if(input$fpkmTransform == 1){
        # log 2 transform. 
        datRes <- mutate(datRes, 
                         value_1 = log2(1+value_1), 
                         value_2 = log2(1+value_2))
      }
      if(input$fpkmTransform == 2){
        # log 10 transform. 
        datRes <- mutate(datRes, 
                         value_1 = log10(1+value_1), 
                         value_2 = log10(1+value_2))
      }
      # no action if fpkmTransform =3
      
      # order genes. 
      if(input$orderGenes == 1){
        datRes <- datRes[(order(datRes$q_value)),]
      }else if(input$orderGenes == 2 & input$pseudoFC){
        datRes <- datRes[(order(abs(datRes$pseudoFC), decreasing=T)),]
      }else{
        datRes <- datRes[(order(abs(datRes$log2.fold_change.), decreasing=T)),]
      }
      
    })
    
    # make nicer table for display. 
    # TO DO : make a submit button for this, or only show when you navigate to data page. 
    
    tableOut <- reactive({
      tabl <- dataInput()[c(3,17,16,18,13,14,2,4)]
      isolate({
        xlab <- paste(input$xAxis, "FPKM")
        ylab <- paste(input$yAxis, "FPKM")
        fclab <- paste("log2 FC", input$yAxis, "rel", input$xAxis)
        
        colnames(tabl)<-c("Gene Symbol", ylab, xlab, fclab, "q value", "significant", "Ensembl ID", "Locus")
        tabl
      })
    })
    
    # making the plots interactive causes a pretty massive slow down. show just static plot first 
    # then give option of showing interactive plot. 
    
    # TO DO - add progress indicator- currently not working great. 
    
    plots <- reactive({
     
      sampleXname <- input$xAxis
        
      sampleYname <- input$yAxis
        
      
      # determine axis transformation. 
      if(input$fpkmTransform == 1){
        fpkmTransformRes <- "log2 FPKM"
      }else if (input$fpkmTransform == 2){
        fpkmTransformRes <- "log10 FPKM"
      }else{
        fpkmTransformRes <- "FPKM"
      }
      
      
      ylabel <- paste(fpkmTransformRes, input$yAxis, sep=" ")
      xlabel <- paste(fpkmTransformRes, input$xAxis, sep=" ")
      
      
      if((dataInput()[1, "sample_1"] == sampleXname) & (dataInput()[1, "sample_2"] == sampleYname)){
        gp <- ggplot(dataInput(), aes(x=value_1, y=value_2))+geom_point(alpha =.5) + theme_classic()+
          xlab(xlabel) + ylab(ylabel) + scale_x_continuous(expand=c(0,0), limits=c(0,12))+
          scale_y_continuous(expand=c(0,0), limits=c(0,12))
      }else{
        gp <- ggplot(dataInput(), aes(x=value_2, y=value_1))+geom_point(alpha =.5) + theme_classic()+
          xlab(xlabel) + ylab(ylabel) + scale_x_continuous(expand=c(0,0), limits=c(0,12))+
          scale_y_continuous(expand=c(0,0), limits=c(0,12))
      }
    })
    
    gv <- reactive({
      # if not show interactive plot, just return empty gv object.  this helps speed up. 
      #dummy plot, 
      if(!input$showInteract){
        gv <- ggvis()
      }else{
        sampleXname <- input$xAxis
          
        sampleYname <- input$yAxis
        
        # determine axis transformation. 
        if(input$fpkmTransform == 1){
          fpkmTransformRes <- "log2 FPKM"
        }else if (input$fpkmTransform == 2){
          fpkmTransformRes <- "log10 FPKM"
        }else{
          fpkmTransformRes <- "FPKM"
        }
        
        ylabel <- paste(fpkmTransformRes, input$yAxis, sep=" ")
        xlabel <- paste(fpkmTransformRes, input$xAxis, sep=" ")
        
        if((dataInput()[1, "sample_1"] == sampleXname) & (dataInput()[1, "sample_2"] == sampleYname)){
          gv <- ggvis(dataInput(), props(x = ~value_1, y= ~value_2))+layer_point(props(fillOpacity := .5))
        }else{
          gv <- ggvis(dataInput(), props(x = ~value_2, y= ~value_1))+layer_point(props(fillOpacity := .5))
        }
        
        # tool tip text. 
        all_values <- function(x) {
          if(is.null(x)) return(NULL)
          format(x$gene)
          # to give values: 
          #paste0(names(x), ": ", format(x), collapse = "<br />")
        }
        
        
        gv <- gv + guide_axis(type="x", title=xlabel) +
          guide_axis(type="y", title=ylabel) + 
          mark_point(props(size := 50, size.hover := 200, 
                           fillOpacity := .1, fillOpacity.hover :=.5,
                           fill.hover := "red",
                           key := ~gene))+tooltip(all_values) 
      }
      
    })
    
    geneExpression<-reactive({
      expressionPlot <- genePlotting(input$geneId, fpkm)
    })
    
    ################### below are output calls. 
    
    output$plot1 <- renderPlot({
      print(plots())
    })  
    
    observe_ggvis(gv, "myplot", session)
    
    output$titleText <- renderText({
      paste(input$yAxis, "vs", input$xAxis
      ) 
    })
    
    output$titleText2 <- renderText({
      paste0(input$yAxis, " vs ", input$xAxis, ". genes passing filters: ", (nrow(dataInput()))
      ) 
    })
    
    
    output$text1 <- renderText({
      print(paste("Number of genes passing filters:", (nrow(dataInput()))
      ))
    })
    
    output$textforDataTable <- renderText({
      print(paste("showing", (nrow(dataInput())), "genes passing filters"))
    })
    
    output$table1 <- renderTable({
      # TO DO : clean up output table widths. not critical now.  
      tableOut()
    }, digits=4, include.rownames = FALSE)
    
    output$downloadData <- downloadHandler(
      filename = c("data.txt"), 
      content = function(file){
        write.table(tableOut(), file, sep="\t", quote=F, row.names=F)
      })
    
    output$geneExpression <- renderPlot({
      print(geneExpression())
      })
    
    
  })