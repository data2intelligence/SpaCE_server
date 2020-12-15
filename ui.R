library(dplyr)
library(shinycssloaders)

navbarPage("SpaCE: Spatial Cellular Estimator", 
           id="nav",
           
           #tabPanel("Home"),
           
           tabPanel("Run",
                    
                    # Sidebar layout with input and output definitions ----
                    sidebarLayout(
                      
                      # Sidebar panel for inputs ----
                      sidebarPanel(
                        
                        # Input: Select a file ----
                        fileInput("file1", "Upload your ST data set (.csv file):",
                                  multiple = TRUE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        
                        a("Download the example file.",href="ST_mel.csv",download=NA, target="_blank"),
                        
                        # Horizontal line ----
                        tags$hr(),
                        
                        # Submit ----
                        actionButton("go", "Submit")
                        
                      ),
                      
                      # Main panel for displaying outputs ----
                      mainPanel(
                        
                          fluidRow(
                            column(6,
                              tags$h5(tags$b("Cell Fraction")),
                              tags$hr(),
                              column(12,
                                #tableOutput("contents") %>% withSpinner(color="#0dc5c1")
                                selectInput('cellType0', 'Cell Type:', 
                                            choices = c(
                                              "Malignant" = "Malignant",
                                              "B cell" = "B cell",
                                              "CAF" = "CAF",
                                              "Dendritic" = "Dendritic",
                                              "Endothelial" = "Endothelial",
                                              "Macrophage" = "Macrophage",
                                              "Macrophage (M1)" = "Macrophage M1",
                                              "Macrophage (M2)" = "Macrophage M2",
                                              "Macrophage (other)" = "Macrophage other",
                                              "Mast" = "Mast",
                                              "NK" = "NK",
                                              "T CD4" = "T CD4",
                                              "T CD4 (Treg)" = "Treg",
                                              "T CD4 (naive)" = "T CD4 naive",
                                              "T CD4 (heler)" = "T CD4 helper",
                                              "T CD4 (other)" = "T CD4 other",
                                              "T CD8" = "T CD8",
                                              "Unidentifiable" = "Unidentifiable"
                                            ))
                                
                              ),
                              plotOutput("plot1") %>% withSpinner(color="#0dc5c1")
                            ),
                            column(6,
                              tags$h5(tags$b("Cell co-localization")),
                              tags$hr(),
                              column(6,
                                selectInput('cellType1', 'Cell Type (x axis):', 
                                          choices = c(
                                            "Malignant" = "Malignant",
                                            "B cell" = "B cell",
                                            "CAF" = "CAF",
                                            "Dendritic" = "Dendritic",
                                            "Endothelial" = "Endothelial",
                                            "Macrophage" = "Macrophage",
                                            "Macrophage (M1)" = "Macrophage M1",
                                            "Macrophage (M2)" = "Macrophage M2",
                                            "Macrophage (other)" = "Macrophage other",
                                            "Mast" = "Mast",
                                            "NK" = "NK",
                                            "T CD4" = "T CD4",
                                            "T CD4 (Treg)" = "Treg",
                                            "T CD4 (naive)" = "T CD4 naive",
                                            "T CD4 (heler)" = "T CD4 helper",
                                            "T CD4 (other)" = "T CD4 other",
                                            "T CD8" = "T CD8",
                                            "Unidentifiable" = "Unidentifiable"
                                          ))
                              ),
                              column(6,
                                selectInput('cellType2', 'Cell Type (y axis):', 
                                          choices = c(
                                            "B cell" = "B cell",
                                            "CAF" = "CAF",
                                            "Dendritic" = "Dendritic",
                                            "Endothelial" = "Endothelial",
                                            "Macrophage" = "Macrophage",
                                            "Macrophage (M1)" = "Macrophage M1",
                                            "Macrophage (M2)" = "Macrophage M2",
                                            "Macrophage (other)" = "Macrophage other",
                                            "Mast" = "Mast",
                                            "NK" = "NK",
                                            "T CD4" = "T CD4",
                                            "T CD4 (Treg)" = "Treg",
                                            "T CD4 (naive)" = "T CD4 naive",
                                            "T CD4 (heler)" = "T CD4 helper",
                                            "T CD4 (other)" = "T CD4 other",
                                            "T CD8" = "T CD8",
                                            "Unidentifiable" = "Unidentifiable"
                                          )),
                              ),
                              plotOutput("plot2") %>% withSpinner(color="#0dc5c1")
                            )
                          )
                        
                      )
                      # Main panel
                      
                    )

                      
           ),
           
           tabPanel("Help"),
           tabPanel("Contact")
           
)