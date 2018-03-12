library(ggplot2)
library(plotly)
library(shiny)
library(DT)
library(edgeR)
library(statmod)
library(GO.db)
library(org.Mm.eg.db)

navbarPage("AAseq DEseq",
           tabPanel("Volcano Plots",
                    sidebarLayout(
                      sidebarPanel(width = 2,
                                   selectInput("compChoose1", label = "Choose Sample 1",
                                               choices = c( "WT_0h", "WT_3h", "WT_6h",
                                                            "KO_0h", "KO_3h", "KO_6h"),
                                               selected = "WT_0h"),
                                   selectInput("compChoose2", label = "Choose Sample 2",
                                               choices = c( "WT_0h", "WT_3h", "WT_6h",
                                                            "KO_0h", "KO_3h", "KO_6h"),
                                               selected = "KO_0h"),
                                   checkboxInput("scaleCheck", "Scale Points to RPKM", F),
                                   selectInput("geneChoose", label = "Gene for DiffTab Barplot",
                                               choices = colnames(readRDS("./Data/AAseqRPKMmat.RData"))),
                                   checkboxInput("returnPath", "Print Pathways Analysis", F),
                                   checkboxInput("pathCheck", "KEGG Pathway (uncheck for GO)", T)
                      ),
                      mainPanel(
                        tabsetPanel(type = "tabs",
                                    tabPanel(
                                      "Volano Plot", plotlyOutput('volcanoOut', height = 600),
                                      plotlyOutput("barPlotNaive"),
                                      verbatimTextOutput("pathPrint")),
                                    tabPanel("Diff Tables", DT::dataTableOutput("diffTable"),
                                             plotlyOutput("barPlotNaive2"))
                        )))),
           tabPanel("Time Course Plot",
                    sidebarLayout(
                      sidebarPanel(width = 2,
                                   selectInput("TCseq", "Transcript", 
                                               choices = colnames(readRDS("./Data/AAseqRPKMmat.RData"))))
                      ,
                      mainPanel(
                        tabsetPanel(type = "tabs",
                                    tabPanel(
                                      "Time Course Plot", plotOutput('TCplot', height = 600))
                        )
                      )
                    )
           ),
           tabPanel("Export RPKM Data",
                    downloadButton('downloadData', 'Download Data'),
                    DT::dataTableOutput('dataDownload'),
                    verbatimTextOutput('dataSelect'),
                    textOutput('selected'))
)
