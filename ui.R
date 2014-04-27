library(ggvis)
shinyUI(fluidPage(
  titlePanel("iMN paper"), 
  # TO DO: make side bar smaller.
  # TO DO: add progress indicator.
  fluidRow(
    column(3, 
           wellPanel(                
             selectInput("yAxis", 
                         label = "choose sample to display on Y axis", 
                         choices = c("Embryonic MNs", "iMNs", "ESC MNs", "iPSC MNs", 
                                     "ESC", "iPSC", "MEF1", "MEF2", "MEF3", "MEF4"), 
                         selected = "iMNs"), 
             
             selectInput("xAxis", 
                         label = "choose sample to display on X axis", 
                         choices = c("Embryonic MNs", "iMNs", "ESC MNs", "iPSC MNs", 
                                     "ESC", "iPSC", "MEF1", "MEF2", "MEF3", "MEF4"), 
                         selected = "Embryonic MNs"),
             
             sliderInput("minFoldChange", 
                         label= "minimum log2 fold change between samples:", 
                         min = .25, max = 5, value = c(1), step = .2),
             br(),
             checkboxInput("sigOnly", 
                           label = "show only significant genes", 
                           value = T), 
             br(),
            
             helpText("Interactive plots allow 
                      identification of genes by mouse over but may be slower to load"),
             br(),
             
             
             
             
             br(),
             br(),
             br(),
             br(),
             helpText("Additional controls include FDR limit, FPKM filtering and transformation"),
             
            
             
             checkboxInput("optional", label= "show optional controls"),
             conditionalPanel(
               condition = "input.optional == true", 
               
               sliderInput("minSig", 
                           label = "maximal False Discovery Rate", 
                           min = 0, max = .2, value = .05, step = .01),
               br(),
               
               radioButtons("orderGenes",
                            label = "Order genes by:", 
                            choices = list("Significance" =1, "|Fold Change|" =2), 
                            selected =1),
               br(),
               br(),
               numericInput("minXfpkm", 
                            label = "minimum FPKM sample X", value=0),
               numericInput("maxXfpkm", 
                            label = "maximum FPKM sample X", value = 5000),
               numericInput("minYfpkm", 
                            label = "minimum FPKM sample Y", value = 0),
               numericInput("maxYfpkm", 
                            label = "maximum FPKM sample Y", value = 5000), 
               
               br(),
               br(),
               radioButtons("fpkmTransform", 
                            label = "Transformation of FPKM", 
                            choices = list("log2" = 1, "log10" = 2 , "no transformation" = 3), 
                            selected = 1), 
               br(),
               checkboxInput("pseudoFC", 
                             label = "Calculate Pseudo fold change", 
                             value = T)
             )
           )
    ),
    column(9, 
           tabsetPanel(
             tabPanel("Plot", 
                      h1(textOutput("titleText")), 
                      h3(textOutput("text1")),
                      plotOutput("plot1", height="500px", width="500px")
                      ),
             tabPanel("Interactive Plot",
                      h3(textOutput("titleText2")), 
                      actionButton("interact", "Refresh Interactive Plot"),
                      br(),
                      ggvis_output(plot_id="myplot")
                      ),
             tabPanel("Data", 
                      # TO DO position on right
                      downloadButton("downloadData", "Save selected Data as CSV"),
                      tableOutput("table1")
             )
           )
    )
  )
))
