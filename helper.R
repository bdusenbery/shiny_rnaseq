genePlotting <- function(GOI, fpkmDf, ...){
  res <- fpkmDf[(fpkmDf$gene %in% GOI),]
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

