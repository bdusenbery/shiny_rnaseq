genePlotting <- function(GOI, fpkmDf, ...){
  res <- fpkmDf[(fpkmDf$ensembl_id %in% GOI),]
  # determine max value
  max <- max(unlist(res$hi_error))
  plot <- ggplot(res, aes(x=sample_good_name, y=fpkm, fill=sample_good_name))+geom_bar(stat="identity")
  plot <- plot+geom_errorbar(aes(min=lo_error, max=hi_error), width=.2)+xlab("")+ylab("FPKM")+
    theme_bw()+scale_y_continuous(expand =c(0,0), limits=c(0, (1.1*max)))+
    scale_fill_manual(values=c("#666699", "#0DB14B", "#008ccf", "#118697", 
                               "#744c28", "#8a5d3b", "#9a8478", "#594a41", 
                               "#d8403f", "#f7931d"))+
    theme(legend.position="none", axis.text.x=element_text(size=12, angle=45, hjust=1))
  return(plot)
}

sigMatrix <- function(GOI, diffDF){
  res <- diffDF[(diffDF$gene_id %in% GOI),]
  res2 <- res
  colnames(res2)[5] <- "sample_2"
  colnames(res2)[6] <- "sample_1"
  res <- rbind(res, res2)
  res$sample_2 <- with(res, factor(sample_2, levels=rev(levels(sample_2))))
  sigMatPlot <- ggplot(res, 
                       aes(x=sample_1, y=sample_2, 
                           fill=cut(q_value, breaks=c(0, .05, .1, .2, .5, 1), 
                                    labels=c("<0.05", "0.05-0.10", "0.10-0.20", "0.20-0.50", ">0.50"), 
                                    include.lowest=T),
                           label=round(q_value, 4)))+geom_tile()+geom_tile(color="black", show_guide=FALSE)+
    scale_fill_manual(values=c("#cb181d", "#fb6a4a", "#fcbba1", "#bababa", "#4d4d4d"), name="FDR")+
    theme_bw()+theme(panel.grid=element_blank())+scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+
    scale_size(guide="none")+xlab("")+ylab("")+theme(axis.text.y=element_text(size=12), 
                                                     axis.text.x=element_text(size=12, angle=45, hjust=1))
  return(sigMatPlot)
}

## process isoforms. 
isoformProcess <- function(GOI, isoFPKM){
  res <- filter(isoFPKM, gene_id == GOI)
  res <- group_by(res, tracking_id) 
  res <- mutate(res, avIso=mean(value))
  res <- arrange(res, desc(avIso))
  res <- group_by(res, variable, add=F)
  res <- mutate(res, sumSample = sum(value))
  return(res)
}



# generate isoform plots
transcriptPlot<- function(isoformProcessRes, type="Bar Plot", normalize=1, ...){
  # set order for isoforms 
  isoformOrder <- unique(as.character(isoformProcessRes$tracking_id))
  isoformProcessRes$tracking_id <- factor(isoformProcessRes$tracking_id, isoformOrder, ordered=T)
  
  # use brewer scales if number of isoforms is less than 8
  isoformNo <- length(isoformOrder)
  # bar chart. 
  if(type == "Bar Plot"){
    if(normalize == 1){
      isoPlot <- ggplot(isoformProcessRes, aes(x=sample_good_name, y=(value/sumSample)*100, 
                                               fill=tracking_id))+geom_bar(stat="identity") + 
        ylab("Percent Isoform Expression") +
        xlab("")
      maxVal <- 97
      
    }else{
      isoPlot <- ggplot(isoformProcessRes, aes(x=sample_good_name, y=(value), 
                                               fill=tracking_id))+geom_bar(stat="identity")+
        ylab("Isoform FPKM") + 
        xlab("")
      maxVal <- max(unlist(isoformProcessRes$sumSample))
    }
    isoPlot <- isoPlot + theme_bw() + scale_y_continuous(expand =c(0,0), limits=c(0, maxVal*1.05))+
      theme(axis.text.x =element_text(size =12, angle=45, hjust=1))
  }else{
    if(normalize == 1){
      isoPlot <- ggplot(isoformProcessRes, aes(x=sumSample/2, 
                                               y=value, fill=tracking_id, width=sumSample))+
        geom_bar(position="fill", stat="identity") +theme_bw()+xlab("FPKM")
    }else{
      isoPlot <- ggplot(isoformProcessRes, aes(x=factor(1), 
                                               y=((value/sumSample)*100), fill=tracking_id))+
        geom_bar(stat="identity",width=1) +theme_bw()+theme(axis.text.y=element_blank(), 
                                                            axis.ticks.y=element_blank())+
        xlab("")
      
    }
    isoPlot <- isoPlot + facet_wrap(facets=(~sample_good_name), ncol=4)+
      coord_polar("y")+ylab("")+
      theme(axis.text.x=element_blank(), panel.grid.major.x=element_line(color="#878787"),
            strip.background = element_blank(), panel.border=element_blank(),
            strip.text=element_text(size=12))
  }
  if(isoformNo <= 8){
    isoPlot <- isoPlot + scale_fill_brewer(type="qual", palette=7, name="Transcript")
  }else{
    isoPlot <- isoPlot +scale_fill_hue(h=c(10, 340), name="Transcript")
  }
  return(isoPlot)
}

methPlot <- function(res){
  methP <- ggplot(res, aes(x=variable, y=value, fill=variable))+geom_bar(stat="identity")
  methP + scale_y_continuous(expand=c(0,0), limits=c(0,1.05))+theme_bw()+facet_grid(position~.)+
    scale_fill_manual(values=c("#666699", "#0DB14B", "#008ccf", "#118697", 
                               "#744c28", "#8a5d3b", "#9a8478", "#594a41", "tan",
                               "#d8403f", "#f7931d"))+xlab("")+ylab("Methylation Ratio")+
    geom_hline(y=.2, color="#bdbdbd", linetype="dashed")+
    geom_hline(y=.6, color="#bdbdbd", linetype="dashed")+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12), legend.position="none")
}