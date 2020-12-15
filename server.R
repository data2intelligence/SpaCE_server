library(BiocManager)
options(repos = BiocManager::repositories())

options(shiny.maxRequestSize = 30*1024^2)

ppHeatmap <- function(visual.matrix,item)
{
  library(ggplot2)
  coordi <- round(t(matrix(as.numeric(unlist(strsplit(colnames(visual.matrix),"x"))),nrow=2)),0)
  fig.df <- data.frame(x=coordi[,1],y=coordi[,2],value=visual.matrix[item,])
  
  p <- ggplot(fig.df,aes(x=x,y=y))+
    geom_point(aes(colour=value),size= 1.8)+
    scale_colour_gradientn(colours = c("blue","yellow","red"))+
    theme_bw()+ 
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid=element_blank(),
      panel.border = element_blank()
    )
  p
}



function(input, output) {
  
  STfile <- eventReactive(input$go, {
    req(input$file1)
  })
  
  res <- reactive({
    
    source("SpaCE.R"); load("Ref.rda")
    ST <- as.matrix(read.csv(STfile()$datapath,row.names=1))
    colnames(ST) <- gsub("X","",colnames(ST))
    
    Signature <- Ref$sigGenes[!names(Ref$sigGenes)%in%"Malignant"]
    malProp <- inferMalignant(ST,Signature,"./infercnv/","x")  # ref, cluster
    
    #malProp <- rep(0,dim(ST)[2])
    #names(malProp) <- colnames(ST)
    
    propMat <- SpatialDeconv(
      ST=ST,
      Ref=Ref,
      malInput=TRUE,
      malProp=malProp,
      Unidentifiable=TRUE, 
      MacrophageOther=TRUE,
      TCD4Other=TRUE
    )
    
    propMat
  })

  
  output$contents <- renderTable(rownames = TRUE,{
    return(res()[,1:3])
  })
  
  output$plot1 <- renderPlot({
    ppHeatmap(res(),input$cellType0)
  })
  
  output$plot2 <- renderPlot({
    plot(t(res()[c(input$cellType1,input$cellType2),]))
  })
  
  
}