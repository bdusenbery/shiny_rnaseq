library(ggvis)
library("shinyIncubator")
shinyUI(navbarPage(
  "iMN paper resource", 
  # About panel -------------------------------------------------------------------
  tabPanel("About", 
           fluidRow(
             column(9,
                    includeMarkdown("about.md")
             ),
             column(3
                    
                    # extra spacing)
             )
           )
  ),
  # comparison Explorer controls -----------------
  tabPanel("Comparison Explorer",
           
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
                      br(),
                      helpText("Additional controls include FDR limit, 
                               FPKM filtering and transformation"),
                      checkboxInput("optional", label= h5("show optional controls")),
                      conditionalPanel(
                        condition = "input.optional == true", 
                        br(),
                        br(),
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
             # comparison explorer outputs --------
             column(9, 
                    tabsetPanel(
                      # comparison explorer static plot -----------
                      tabPanel("Plot", 
                               h1(textOutput("titleText")), 
                               h3(textOutput("text1")),
                               plotOutput("plot1", height="500px", width="500px")
                      ),
                      # comparison explorer interactive plot -----------
                      tabPanel("Interactive Plot", 
                               fluidRow(
                                 column(4, 
                                        checkboxInput("showInteract", "show interactive plot", value=FALSE)
                                 ),
                                 column(8, 
                                        # use an action button here to isolate the interactive plot
                                        # which requires significant loading/comput time. 
                                        helpText("Interactive plot function allows individual genes to 
                                                 identified by mouse over, but reduces app speed overall.
                                                 Uncheck box to improve performance.  It may take a moment 
                                                 for the interactive plot to appear.")
                                 ),
                                 br(),
                                 fluidRow(
                                   textOutput("interactText"),
                                   #plotOutput("interactivePlot")
                                   ggvis_output(plot_id="myplot")
                                 )
                               )
                      ), 
                      # comparison exploerer data download ----------
                      tabPanel("data",
                               fluidRow(
                                 column(8, 
                                        textOutput("textforDataTable")), 
                                 column(4,
                                        downloadButton("downloadData", "Save selected Data as CSV")
                                 )
                               ),
                               br(),
                               tableOutput("table1")
                      )
                    )
             )
           )
  ),
  # gene Explorer -----------
  tabPanel("Gene Explorer",
           fluidRow(
             column(3, 
                    wellPanel(
                      helpText("Eventually will hopefully allow autocomplete of gene names, 
                               but for now you must enter a proper ID, case sensitive"),
                      # TO DO : better gene input see shiny sky or adapt selectize js. 
                      # TO DO : error handling for no gene and duplicate genes.
                      # TO DO : provide ensemble ID input 
                      # TO DO : provide submit button. 
                      # TO DO : consider making this float.
                      textInput("geneId", label="Select gene", value="Mnx1")
                      
                    )), 
             # gene Explorer expression plots ------------
             column(9, 
                    tabsetPanel(
                      tabPanel("Gene Expression", 
                               # TO DO : put gene id in a box, make ensembl linkout
                               wellPanel(
                                 h3(textOutput("geneInfo1"))
                               ),
                               h4("Gene expression across samples"),
                               plotOutput("geneExpression"),
                               br(),
                               h4("Significance Matrix between samples"),
                               plotOutput("sigMat",  height="500px", width="600px")
                               
                      ),
                      # Gene Explorer isoform plots --------------
                      tabPanel("Isoform Expression", 
                               # TO DO : consider printing out isoform expression or providing table- 
                               # currently can't copy and paste isoform ids. 
                               # alternatively it would be awesome to show structure and even more awesome
                               # if we could do this in the same colors as the pies 
                               wellPanel(
                                 fluidRow(
                                   column(3, 
                                          # select type of plot to show for isoform plots
                                          selectInput(inputId="showIsoType", 
                                                      label="Show Isoforms as:", 
                                                      choices= c("Bar Plot", "Pie Chart"), 
                                                      selected="Pie Chart")
                                   ),
                                   column(3, 
                                          # select if isoforms should be normalized to total expression
                                          radioButtons(inputId="normalizeIsoforms", 
                                                       label="Normalize Isoforms to total gene expression", 
                                                       choices=list("yes" = 1, "no" = 2), 
                                                       selected = 1)
                                   ),
                                   column(3, 
                                          selectInput(inputId="numberIso", 
                                                      label="Maximum number of most abundant isoforms to show:", 
                                                      choices = c(1:20), 
                                                      selected = 8)
                                   ),
                                   column(3, 
                                          helpText("Normalizing to total gene expression shows relative levels
                                                 of isoforms, otherwise percents of each isoform is shown. 
                                                 Most abundant isoforms are determined as mean across all samples")
                                   )
                                 )
                               ),
                               fluidRow(
                                 h3(textOutput("isoformReportingText")),
                                 br(),
                                 plotOutput("isoformPlot", height="500px", width="600px")
                                 
                               )
                      ),
                      # Gene explorer promoter methylation. 
                      tabPanel("Promoter Methylation", 
                               fluidRow(
                                 wellPanel(
                                   h3(textOutput("methylGeneInfo")),
                                   h4(textOutput("methylGenePosition"))
                                 )),
                               fluidRow(
                                 plotOutput("methylationPlot", height="500px", width="600px")
                               )
                      )
                    ))
             
           ))
))



