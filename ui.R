#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

fluidPage(
    titlePanel("Visualize Spatial Transcriptomics Data (Private Website of Jiang's Lab)"),
    
    fluidRow(
        column(3,wellPanel(
            selectInput('st', 'Spatial Transcriptomics (ST) Data Set:', sts, selected = sts[3]),
            selectInput('sampleName', 'Sample from the selected ST Data Set:', 'start sample'),
            
            tags$h5(tags$b("Parameter for module 1 (Cell)")),
            wellPanel(
                selectInput('sc', 'Single Cell (SC) Data Set:', scs, selected = scs[1]),
                #radioButtons("group", "Signature Matrix",choices = c("Cell Type" = "uniformCellType","Cell Type (Sub)" = "uniformCellTypeSub"), inline = TRUE),
                selectInput('cellType', 'Cell Type:', 
                            choices = c(
                                "B cell" = "B cell",
                                "CAF" = "CAF",
                                "Dendritic" = "Dendritic",
                                "Endothelial" = "Endothelial",
                                "Macrophage" = "Macrophage",
                                "Macrophage (M1)" = "M1",
                                "Macrophage (M2)" = "M2",
                                "Macrophage (other)" = "Mother",
                                "Mast" = "Mast",
                                "NK" = "NK",
                                "T CD4" = "T CD4",
                                "T CD4 (Th)" = "Th",
                                "T CD4 (Treg)" = "Treg",
                                "T CD4 (other)" = "Tother",
                                "T CD8" = "T CD8",
                                "unknownCell" = "unknownCell"
                            )),
                tags$style(type='text/css', '#txtout1 {color: red;}'),
                verbatimTextOutput("txtout1"),
                radioButtons("colorScale", "Color",choices = c("Green" = "color1","Red and Blue" = "color2"), inline = TRUE)
            ),
            
            tags$h5(tags$b("Parameter for module 2 (Gene)")),
            wellPanel(
                textInput("posGenes", "Positive Genes (separated by space):", "MS4A1"),
                textInput("negGenes", "Negative Genes (separated by space):", ""),
                verbatimTextOutput("txtout2"),
                tags$style(type='text/css', '#txtout3 {color: red;}'),
                verbatimTextOutput("txtout3")
            ),
            
            tags$h5(tags$b("Parameter for module 3 (Clustering)")),
            wellPanel(
                sliderInput("km", "K-means", min=2, max=10, value=4, step=1)
            )
        )),
        
        column(9,wellPanel(
            fluidRow(
                column(6,wellPanel(
                    h4("Module 1: Cell Type Abundance"),
                    plotOutput("plot1",height = figHeight, width = figWidth)
                )),
                column(6,wellPanel(
                    h4("Summary for module 1"),
                    plotOutput("plot1.1",height = figHeight, width = figWidth)
                ))
            ),
            fluidRow(
                column(6,wellPanel(
                    h4("Module 2: Gene Expression"),
                    plotOutput("plot2",height = figHeight, width = figWidth)
                )),
                column(6,wellPanel(
                    h4("Module 1 (Cell) vs Module 2 (Gene)"),
                    plotOutput("plot4",height = figHeight, width = figWidth)
                ))
            ),
            fluidRow(
                column(6,wellPanel(
                    h4("Module 3: Clustering"),
                    plotOutput("plot3",height = figHeight, width = figWidth)
                )),
                column(6,wellPanel(
                    h4("Module 1 (Cell) vs Module 3 (Clustering)"),
                    plotOutput("plot5",height = figHeight, width = figWidth)
                    #h4("Module 2 (Gene) vs Module 3 (Clustering)"),
                    #plotOutput("plot6",height = figHeight, width = figWidth)
                ))
            ),
        ))
        
    )
)
