#
# This is a Shiny web application called T-rex. You can run the application by clicking
# the 'Run App' button above.
#
# Authors: David Bergenholm, Christoph Börlin, Petter Holland & Jens Nielsen
# Chalmers University of Technology
# Department of Biology and Bioengineering

library(rsconnect)




suppressMessages(library(shiny))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(DT))
suppressMessages(library(gdata))
suppressMessages(library(gplots))
suppressMessages(library(gtools))
suppressMessages(library(scales))
suppressMessages(library(GGally))
suppressMessages(library(network))
suppressMessages(library(MASS))
suppressMessages(library(ggpubr))
suppressMessages(library(cluster))
suppressMessages(library(factoextra))
suppressMessages(library(biomaRt))
suppressMessages(library(shinyWidgets))
suppressMessages(library(Biostrings))
suppressMessages(library(IRanges))
suppressMessages(library(shinycssloaders))
# suppressMessages(library(NbClust))
suppressMessages(library(sna))


# Define UI for application that draws a histogram
fluidPage(
  
  # Application title
  tabsetPanel(type="tabs",
              
              tabPanel("Welcome",
                       fluidRow(
                         column(width=12, align="center", 
                                h1(" "),
                                h3("Welcome to"),
                                h2("T-rEx"),
                                h1(" "),
                                h4(tags$em("S. cerevisiae")," Transcription factor Explorer"),
                                h6("David Bergenholm, Christoph Börlin, Petter Holland, & Jens Nielsen"),
                                imageOutput("sctrex",height=120) %>% withSpinner(color="#0dc5c1")
                  
                         ))),
              tabPanel("Transcription Factor Summary",
                       sidebarLayout(
                         sidebarPanel(width=2,
                           fluidRow("Select Transcription Factor",
                                   column(width=6,
                                          radioButtons("TFs", "",
                                                 list("Cat8"='Cat8', 
                                                      "Cbf1"='Cbf1', 
                                                      "Ert1"='Ert1', 
                                                      "Gcn4"='Gcn4',
                                                      "Gcr1"='Gcr1', 
                                                      "Gcr2"='Gcr2', 
                                                      "Hap1"='Hap1', 
                                                      "Ino2"='Ino2', 
                                                      "Ino4"='Ino4', 
                                                      "Leu3"='Leu3',
                                                      "Oaf1"='Oaf1', 
                                                      "Pip2"='Pip2', 
                                                      "Rds2"='Rds2',	
                                                      "Rgt1"='Rgt1',	
                                                      "Rtg1"='Rtg1',	
                                                      "Rtg3"='Rtg3',	
                                                      "Sip4"='Sip4',	
                                                      "Stb5"='Stb5',	
                                                      "Sut1"='Sut1',	
                                                      "Tye7"='Tye7')))
                                           ,
                                    column(width=12,
                                           selectInput("dataseries", "Select Data",c("GEM: Peaks"), selected="GEM: Peaks", selectize=TRUE),# Could include sum of reads but will not be included in this version
                                    h4("Dowload Data"),
                                    h5("Peak list contains the sum of all peaks identied at a given promoter for each condition."),
                                    h5("Peak list advanced contains the gene, chromosome, strand, peak position, strength of binding and the condition"),
                                    downloadButton("downloadPeakData", "Download Peaklist"),
                                    downloadButton("downloadPeakDataAd", "Download Peaklist advanced"))
                                   
                                    
                           )),
                         mainPanel(
                           fluidRow( 
                                           
                                           h4("Targets"),
                                           h6(dataTableOutput("table1",height=200 ))),
                           fluidRow(
                             
                                    
                                    h4("Consensus Motif"),
                                    h6("The consensus motif is genreated for each condition "),
                                    column(width=3,
                                           h5("Glu-lim"),
                                           imageOutput("motif_glu",height=120)),
                                    column(width=3,
                                           h5("Nit-lim"),
                                           imageOutput("motif_nit",height=120)),
                                    column(width=3,
                                           h5("Eth-lim"),
                                           imageOutput("motif_eth",height=120)),
                                    column(width=3,
                                           h5("Ana-lim"),
                                           imageOutput("motif_ana",height=120)),
                                    h6(""),
                                    column(width=3,
                                           
                                           imageOutput("motif_glu_seq",height=120)),
                                    column(width=3,
                                           
                                           imageOutput("motif_nit_seq",height=120)),
                                    column(width=3,
                                           
                                           imageOutput("motif_eth_seq",height=120)),
                                    column(width=3,
                                           
                                           imageOutput("motif_ana_seq",height=120)),
                                    
                                    h4("Peak Distribution Profile"),
                                    h6("The peak distribution profile is genreated for each condition "),
                                    column(width=3,
                                           h5("Glu-lim"),
                                           imageOutput("peakdist_glu",height=200)),
                                    column(width=3,
                                           h5("Nit-lim"),
                                           imageOutput("peakdist_nit",height=200)),
                                    column(width=3,
                                           h5("Eth-lim"),
                                           imageOutput("peakdist_eth",height=200)),
                                    column(width=3,
                                           h5("Ana-lim"),
                                           imageOutput("peakdist_ana",height=200))
                                    
                            )
                         )
                       )
              ),
              tabPanel("Transcription Factor Binding Data",
                       fluidRow(
                         column(width=2,class="well",

                                textInput("text", label = h4("Select Gene"), value = "ACC1"),
                                verbatimTextOutput("GeneInfo"),
                                h4("Select Transcription Factor"),
                                column(width=6,
                                checkboxInput("Cat8","Cat8", value = F),
                                checkboxInput("Cbf1","Cbf1", value = F),
                                checkboxInput("Ert1","Ert1", value = F),
                                checkboxInput("Gcn4","Gcn4", value = F),
                                checkboxInput("Gcr1","Gcr1", value = F),
                                checkboxInput("Gcr2","Gcr2", value = F),
                                checkboxInput("Hap1","Hap1", value = F),
                                checkboxInput("Ino2","Ino2", value = F),
                                checkboxInput("Ino4","Ino4", value = F),
                                checkboxInput("Leu3","Leu3", value = F),
                                checkboxInput("Oaf1","Oaf1", value = F)),
                                column(width=6,
                                checkboxInput("Pip2","Pip2", value = F),
                                checkboxInput("Rds2","Rds2", value = F),
                                checkboxInput("Rgt1","Rgt1", value = F),
                                checkboxInput("Rtg1","Rtg1", value = F),
                                checkboxInput("Rtg3","Rtg3", value = F),
                                checkboxInput("Sip4","Sip4", value = F),
                                checkboxInput("Stb5","Stb5", value = F),
                                checkboxInput("Sut1","Sut1", value = F),
                                checkboxInput("Tye7","Tye7", value = F)),
                                actionButton("Load", "Load Data"),
                                
                                column(width=12,
                                textInput("motifstr", label = h4("Motif"), value = ""),
                                checkboxInput("motiffinder","Search", value = F),
                                verbatimTextOutput("info"),
                                h4("Sequence"),
                                textInput("SeqMin", label = h5("Sequence From"), value = ""),
                                textInput("SeqMax", label=h5("Sequence To"),value = ""),
                                textOutput("Sequnce_out"),
                                checkboxInput("SeqFind","Search", value = F),
                                selectInput("ReadsData", "Download Data",c("Glu","Nit","Eth","Ana"), selected="Glu", selectize=TRUE),
                                h5("Download the displayed reads"),
                                downloadButton("downloadreadsData", "Download Reads"))
                                
                                
                         ),   
                         
                         
                         column(width=10, 
                                
                                column(width=9, 
                                       h4("Transcription factor binding profile"),
                                       h5("Glu-lim"),
                                       plotOutput("plot1", 
                                                  height = 180, 
                                                  dblclick = "plot1_dblclick", 
                                                  hover ="plot1_hover",
                                                  brush = brushOpts(id = "plot1_brush",  
                                                                    resetOnNew = TRUE)) %>% withSpinner(color="#0dc5c1"),
                                       
                                       h5("Nit-lim"),
                                       plotOutput("plot2", 
                                                  height = 180, 
                                                  dblclick = "plot1_dblclick", 
                                                  hover ="plot1_hover",
                                                  brush = brushOpts(id = "plot1_brush",  
                                                                    resetOnNew = TRUE))%>% withSpinner(color="#0dc5c1"),
                                       
                                       h5("Eth-lim"),
                                       plotOutput("plot3", 
                                                  height = 180, 
                                                  dblclick = "plot1_dblclick", 
                                                  hover ="plot1_hover",
                                                  brush = brushOpts(id = "plot1_brush",  
                                                                    resetOnNew = TRUE))%>% withSpinner(color="#0dc5c1"),
                                       h5("Ana-lim"),
                                       plotOutput("plot4", 
                                                  height = 180, 
                                                  dblclick = "plot1_dblclick", 
                                                  hover ="plot1_hover",
                                                  brush = brushOpts(id = "plot1_brush",  
                                                                    resetOnNew = TRUE))%>% withSpinner(color="#0dc5c1")),
                                
                                column(width=1, 
                                       h4("Legend"),
                                       checkboxInput("TF_BS","TF BS", value = F),
                                       checkboxInput("GTF","TATA", value = F),
                                       checkboxInput("yranges","Fixed y-axis", value = F),
                                       plotOutput("plot5",
                                                  height = 800,
                                                  dblclick = "plot1_dblclick",
                                                  hover ="plot1_hover",
                                                  brush = brushOpts(id = "plot1_brush",
                                                                    resetOnNew = TRUE))))
                       )
                       #)
                       #)
              ),
              
              tabPanel("Transcription Factor Analysis",
                       sidebarLayout(
                         sidebarPanel(
                                      h3("Dataset"),
                                      h5("The analysis is based on peak identification"),
                                      selectInput("val1", "Exclude or Include data",c("Exclude dubious and non-Yeast 8","Include dubious and non-Yeast 8"), selected="Exclude dubious and non-Yeast 8", selectize=TRUE),
                                      textInput("goterm", "GO-term", value="Amino acid"),
                                      h6("To include more GO-terms add + inbetween."),
                                      selectInput("test", "Test",c("Fisher","Heatmap","Network","Cluster","Linear Model"), selected="Fisher", selectize=TRUE),
                                      textOutput("selectInput"),
                                      conditionalPanel(
                                        condition = "input.test == 'Cluster'",
                                        sliderInput("slider1", label = h5("Number of cluster"), min = 1, 
                                                    max = 10, value = 5)),
                                      
                                      actionButton("search", "Search"),
                                      h4("The selected GO-terms from search result"),
                                      dataTableOutput("goterms"),
                                      selectInput("StatData", "Download Data",c("Glu","Nit","Eth","Ana"), selected="Glu", selectize=TRUE),
                                      downloadButton("downloadStatData", "Download")
                         ), 
                         
                         mainPanel(
                           fluidRow(
                             column(width = 12, 
                                    column(width=6,
                                           h4("Glu-lim"),
                                           plotOutput("plot1_p3", 
                                                      height = 350  
                                           )%>% withSpinner(color="#0dc5c1")),
                                    
                                    column(width = 6, 
                                           h4("Nit-lim"),
                                           plotOutput("plot2_p3", 
                                                      height = 350 
                                           )%>% withSpinner(color="#0dc5c1")),
                                    
                                    
                                    column(width=6,
                                           h4("Eth-lim"),
                                           plotOutput("plot3_p3", 
                                                      height = 350 
                                           )%>% withSpinner(color="#0dc5c1")),
                                    
                                    column(width = 6, 
                                           h4("Ana-lim"),
                                           plotOutput("plot4_p3", 
                                                      height = 350 
                                           )%>% withSpinner(color="#0dc5c1"))))
                         )
                       )
              )
  )
)
