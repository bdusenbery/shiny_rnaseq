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
isoFpkm <- readRDS("data/isoform_fpkm.rds")
methDat <- readRDS("data/methylation.rds")
geneIds <- readRDS("data/mouseGeneInfo.rds")


source("helper.R")


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
    
    # gene id input. if input is gene symbol, get ensembl id. 
    # if there are multiple ensembl ids, warn. 
    geneData <- reactive({
      # depend on get gene button. 
      input$update
      isolate({
      if(input$idType == "Official Gene Symbol"){
        geneDat <- geneIds[(geneIds$mgi_symbol %in% input$geneId),]
      }else{
        geneDat <- geneIds[(geneIds$ensembl_gene_id %in% input$geneId),]
      }
      })
    })
    
    # if geneData returns a single entry, show gene information,
    # if more than one, provide info. 
    # if no matches, suggest. 
    # TO DO: make formating a little nicer. 
    geneInfo <- reactive({
      numberGenesFound <- nrow(geneData())
      if (numberGenesFound == 0){
        geneInfo <- "No gene found, did you mean:"
      } else if (numberGenesFound == 1){
        geneInfo <- paste(geneData()$mgi_symbol, 
                           "--", geneData()$description, 
                          "ensembl link:", geneData()$ensembl_gene_id,
                           "entrez gene link:", geneData()$entrezgene 
                           )
      }else{
        dupids <- paste(c(as.character(geneData()$ensembl_gene_id)), collapse=", ")
        
        geneInfo <- paste("Multiple genes found with gene symbol:", unique(geneData()$mgi_symbol), 
                          "Select one Ensembl Id:", dupids)
      }
       })

    # gene expression bar plot 
    geneExpression<-reactive({
      geneId <- geneData()$ensembl_gene_id
      if(length(geneId) == 1){
      expressionPlot <- genePlotting(geneId, fpkm)
      }else{
        expressionPlot <- qplot(1:10, 1:10, geom="blank")+
          geom_text(aes(x=5, y=5), label="Select a new\ngene ID", 
                    size=14, color="#bdbdbd")+theme_classic()+
          theme(axis.text=element_blank(), axis.ticks= element_blank(), 
                line=element_blank())+xlab("")+ylab("")
      }
    })
    
    # gene expression significance matrix
    sigMatRes <- reactive({
      geneId <- geneData()$ensembl_gene_id
      if(length(geneId) == 1){
      sigMatrixRes <- sigMatrix(geneId, dat)
      }else{
        expressionPlot <- qplot(1:10, 1:10, geom="blank")+
          geom_text(aes(x=5, y=5), label="Select a new\ngene ID", 
                    size=14, color="#bdbdbd")+theme_classic()+
          theme(axis.text=element_blank(), axis.ticks= element_blank(), 
                line=element_blank())+xlab("")+ylab("")
      }
    })
    
    ###### controls for isoform expression. 
    processedIso <- reactive({
      geneId <- geneData()$ensembl_gene_id
      if(length(geneId) == 1){
      processedIso <- isoformProcess(geneId, isoFpkm)
      }
    })
    
   
    
    # return total number of isoforms, 
    isoNumbText <- reactive({
      numIso <- nrow(processedIso())/10
      if(numIso <= input$numberIso){
        isotext <- paste("Showing", numIso, "of", numIso, "observed isoforms")
      }else{
        isotext <- paste("Showing", input$numberIso, "of", numIso, "observed isoforms")
      }
      })
    
    # if total number of isoforms exceeds selected input, limit to selected. 
    processedIsoFinal <- reactive({
      numIso <- nrow(processedIso())/10
      if(numIso <= input$numberIso){
        final <- processedIso()
      }else{
        # rows to keep is * 10 because there are 10 samples per isoform. 
        rowsKeep <- input$numberIso * 10
        final <- processedIso()[1:rowsKeep]
      }
    })
    
    isoformTableDat <- reactive({
      show <- processedIsoFinal()[c("tracking_id", "length")]
      show <- show[!duplicated(show$tracking_id),]
      show$transcript <- show$tracking_id
      show <- show[c("transcript", "length")]
    })
    
    # generate isoform plot. 
    isoPlot <- reactive({
      geneId <- geneData()$ensembl_gene_id
      if(length(geneId) == 1){
      isoPlot <- transcriptPlot(processedIsoFinal(), input$showIsoType, input$normalizeIsoforms)
      }else{
        isoPlot <- qplot(1:10, 1:10, geom="blank")+
          geom_text(aes(x=5, y=5), label="Select a new\ngene ID", 
                    size=14, color="#bdbdbd")+theme_classic()+
          theme(axis.text=element_blank(), axis.ticks= element_blank(), 
                line=element_blank())+xlab("")+ylab("")
      }
    })
    
    # select methylation data. 
    # ensembl ids in methylation data are more difficult. 
    # also, sometimes there are multiple promoter regions associated with one gene: 
    methData <- reactive({
      geneId <- geneData()$mgi_symbol
      
        res <- methDat[(methDat$geneNames %in% geneId),]
      
    })
    
    #output text for methylation data gene name. 
    
    #output text for methylation gene postion. 
    methtext2 <- reactive({
      if(nrow(methData()) > 0){
      paste("Promoter location (mm9 assembly) :", as.character(unique(methData()$position)))
      }
      })
    #output plot for methylation data. 
    methPlotOut <- reactive({
      if(nrow(methData() > 0)){
        methPlotOut <- methPlot(methData()) 
      }else{
        methPlotOut <- qplot(1:10, 1:10, geom="blank")+
          geom_text(aes(x=5, y=5), label="No methylation data\nfor this gene", 
                    size=14, color="#bdbdbd")+theme_classic()+
          theme(axis.text=element_blank(), axis.ticks= element_blank(), 
                line=element_blank())+xlab("")+ylab("")
      }
      
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
   
    # Gene Info for gene expression page:
    output$geneInfo1 <- renderText({
      geneInfo()
    })
    
    # gene expression bar plot
    output$geneExpression <- renderPlot({
      print(geneExpression())
      })
    
    # gene sig matrix
    output$sigMat <- renderPlot({
      print(sigMatRes())
    })
    
    # gene data for isoform page. 
    output$isoformStatus <- renderText({
      geneInfo()
    })
    
    # isoform data for table. 
    output$isoformTable <- renderTable({
      isoformTableDat()
    }, include.rownames = FALSE)
    
    # isoform outputs. 
    output$isoformReportingText <- renderText({
      isoNumbText()
    })
    
    # isoform plot 
    output$isoformPlot <- renderPlot({
      print(isoPlot())
    })
    
   
    # gene info display for methylation page. 
    output$methylGeneInfo <- renderText({
      geneInfo()
    })
    
    # methylation text2
    output$methylGenePosition <- renderText({
      methtext2()
    })
    
    # methylation plot
    output$methylationPlot <- renderPlot({
      print(methPlotOut())
    })
    
  })