#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

ppHeatmap <- function(visual.matrix,item,title,legendName,colorScale,imagePath,st,sampleName)
{
    coordi <- round(t(matrix(as.numeric(unlist(strsplit(colnames(visual.matrix),"x"))),nrow=2)),0)
    fig.df_pre <- data.frame(x=coordi[,1],y=coordi[,2])
    fig.df <- rotateDfr(fig.df_pre, degrees[[st]][[sampleName]])
    
    if(st%in%sts.y.reverse) fig.df[,2] <- -fig.df[,2]
    
    r <- png::readPNG(paste0(imagePath,st,"/",sampleName,".png"))
    rg <- grid::rasterGrob(r, width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
    
    fig.df[["Value"]] = visual.matrix[item,]
    
    p <- ggplot(fig.df,aes(x=x,y=y))+ 
        annotation_custom(rg) + # add background image
        xlim(min(fig.df[,1])-margins[[st]][[sampleName]][1],max(fig.df[,1])+margins[[st]][[sampleName]][2])+
        ylim(min(fig.df[,2])-margins[[st]][[sampleName]][3],max(fig.df[,2])+margins[[st]][[sampleName]][4])+
        ggtitle(title)+
        theme_bw()+ 
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid=element_blank(),
            panel.border = element_blank()
        )
    
    if(colorScale=="color1")
    {
        p+geom_point(aes(alpha=Value),size= pointSize[st], color="green")+
        scale_alpha_continuous(name=paste0("                 \n",legendName),range = c(0,1))
    }else if(colorScale=="color2"){
        p+geom_point(aes(colour=Value),size= pointSize[st])+
        scale_colour_gradientn(name=paste0("                 \n",legendName), colours = c("blue","white","red"))
    }else{
        p+geom_point(aes(colour=Value),size= pointSize[st])+
        labs(colour = paste0("                 \n",legendName))
    }
    
}

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    observe({
        sampleNames <- readLines(paste0(namePath,"sampleName_",input$st,".txt"))
        updateSelectInput(session, "sampleName", choices = sampleNames, selected = sampleNames[1])
        
        #cellTypeNames <- readLines(paste0("ST_data/cellType/cellTypeName_",input$group,".txt"))
        #updateSelectInput(session, "cellType", choices = cellTypeNames)
    })
    
    cellMTX <- reactive({
        st <- input$st
        sampleName <- input$sampleName
        sc <- input$sc
        group <- input$group
        
        filePath <- paste0(deconvResPath,st,"/",sampleName,"_Ref_SC","_",sc,".csv")
        
        if(file.exists(filePath))
        {
            visual.matrix <- as.matrix(read.csv(filePath,as.is=T,row.names=1,header=T))
            colnames(visual.matrix ) <- gsub("X","",colnames(visual.matrix))
            visual.matrix <- filterST(st,sampleName,visual.matrix)
            
            list(filePath=filePath, visual.matrix=visual.matrix)
        }
    })
    
    output$txtout1 <- renderText({
        sc <- input$sc
        cellType <- input$cellType
        
        if(is.list(cellMTX()))
        {
            if(!cellType%in%rownames(cellMTX()$visual.matrix))
            {
                paste0("Warning: \nThe SC data set ", sc, " \ndoes not include cell type ", cellType, ".\n")
            }
        }
    })
    
    output$plot1 <- renderPlot({
        st <- input$st
        sampleName <- input$sampleName
        sc <- input$sc
        group <- input$group
        cellType <- input$cellType
        colorScale <- input$colorScale

        if(is.list(cellMTX()) & cellType%in%rownames(cellMTX()$visual.matrix))
        {
            ppHeatmap(cellMTX()$visual.matrix,cellType,cellType,"Prop",colorScale,imagePath,st,sampleName)
        } 
    })
    
    output$plot1.1 <- renderPlot({
        st <- input$st
        sampleName <- input$sampleName
        sc <- input$sc
        group <- input$group
        cellType <- input$cellType
        colorScale <- input$colorScale
        
        if(is.list(cellMTX()) & cellType%in%rownames(cellMTX()$visual.matrix))
        {
            yValue <- sort(cellMTX()$visual.matrix[cellType,], decreasing = T)
            xValue <- 1:length(yValue)
            
            fig.df <- data.frame(x=xValue, y=yValue)
            ggplot(fig.df,aes(x,y))+
                geom_point()+
                xlab("Spots")+
                ylab("Prop")+
                ggtitle(cellType)+
                theme_bw()+ 
                theme(
                    plot.title = element_text(hjust = 0.5)
                )
        } 
    })
    
    exprMTX <- reactive({
        st <- input$st
        sampleName <- input$sampleName
        
        filePath <- paste0(preprocessDataLogCPMPath,st,"/",sampleName,"_logCPM.rda")
        
        if(file.exists(filePath))
        {
            load(filePath)
            visual.matrix <- filterST(st,sampleName,st.matrix.data)
            
            list(filePath=filePath, visual.matrix=visual.matrix)
        }
    })
    
    output$txtout2 <- renderText({
        posGenes <- normalizeInput(input$posGenes)
        negGenes <- normalizeInput(input$negGenes)
        
        paste0(
            "You input: \n",
            length(posGenes)," positive gene(s) \n",
            length(negGenes)," negative gene(s) \n"
        )
    })
    
    output$txtout3 <- renderText({
        posGenes <- normalizeInput(input$posGenes)
        negGenes <- normalizeInput(input$negGenes)
        Genes <- unique(c(posGenes, negGenes))
        
        if(is.list(exprMTX()))
        {
            GenesFilter <- Genes[!Genes%in%rownames(exprMTX()$visual.matrix)]
            
            if(length(GenesFilter)!=0)
            {
                paste0(
                    "Warning:","\n",
                    "The genes (","i.e., ", paste(GenesFilter,collapse =","), ")","\n",
                    "have been excluded because they are","\n",
                    "1) not official symbols OR","\n",
                    "2) not in the expression matrix.","\n"
                )
            }
        }
    })
    
    averMTX <- reactive({
        posGenes <- normalizeInput(input$posGenes)
        negGenes <- normalizeInput(input$negGenes)
        
        if(is.list(exprMTX()))
        {
            posGenes <- posGenes[posGenes%in%rownames(exprMTX()$visual.matrix)]
            negGenes <- negGenes[negGenes%in%rownames(exprMTX()$visual.matrix)]
            
            posGenesL <- length(posGenes)
            negGenesL <- length(negGenes)
            
            averMTX <- exprMTX()$visual.matrix[c(posGenes,negGenes),,drop=F]
            
            if(posGenesL!=0 | negGenesL!=0)
            {
                if(posGenesL!=0&negGenesL==0){
                    temp = colMeans(averMTX[posGenes,,drop=F])
                }else if(posGenesL==0&negGenesL!=0){
                    temp = -colMeans(averMTX[negGenes,,drop=F])
                }else{
                    temp = colMeans(averMTX[posGenes,,drop=F])-colMeans(averMTX[negGenes,,drop=F])
                }
                
                averMTX <- rbind(averMTX,Average=temp)
            }
            
            list(
                posGenes=paste(posGenes,collapse =", "),
                negGenes=paste(negGenes,collapse =", "),
                visual.matrix = averMTX
            )
        }
    })
    
    output$plot2 <- renderPlot({
        st <- input$st
        sampleName <- input$sampleName
        
        if(is.list(averMTX()) & "Average"%in%rownames(averMTX()$visual.matrix))
        {
            title <- paste0("Average (",averMTX()$posGenes,") - Average (",averMTX()$negGenes,")")
            ppHeatmap(averMTX()$visual.matrix,"Average",title,"Log2CPM","color2",imagePath,st,sampleName)
        }
    }) 
   
    clusterMTX <- reactive({
        st <- input$st
        sampleName <- input$sampleName
        
        filePath <- paste0(clusteringPath,st,"/",sampleName,".csv")
        
        if(file.exists(filePath))
        {
            visual.matrix <- as.matrix(read.csv(filePath,as.is=T,row.names=1,header=T))
            colnames(visual.matrix) <- gsub("X","",colnames(visual.matrix))
            visual.matrix <- filterST(st,sampleName,visual.matrix)
            visual.matrix <- apply(visual.matrix, c(1,2), as.character)
            
            list(filePath=filePath, visual.matrix=visual.matrix)
        }
    })
    
    output$plot3 <- renderPlot({
        st <- input$st
        sampleName <- input$sampleName
        km <- as.character(input$km)

        if(is.list(clusterMTX()))
        {
            ppHeatmap(clusterMTX()$visual.matrix, km, paste0("km = ",km), "Cluster", "cluster", imagePath,st,sampleName)
        }
    })
    

    output$plot4 <- renderPlot({
        st <- input$st
        sampleName <- input$sampleName
        sc <- input$sc
        group <- input$group
        cellType <- input$cellType
        km <- as.character(input$km)
        
        if(is.list(cellMTX()) & is.list(averMTX()) & is.list(clusterMTX()) & 
           cellType%in%rownames(cellMTX()$visual.matrix) & "Average"%in%rownames(averMTX()$visual.matrix))
        {
            olp <- intersect(colnames(cellMTX()$visual.matrix),colnames(averMTX()$visual.matrix))

            xValue <- averMTX()$visual.matrix["Average",olp]
            yValue <- cellMTX()$visual.matrix[cellType,olp]
            zValue <-clusterMTX()$visual.matrix[km,olp]
            
            ccdf <- data.frame(xValue=xValue,yValue=yValue,zValue=zValue)
            
            suppressWarnings(cor_obj <- cor.test(xValue,yValue,method = "spearman"))
            rv <- signif(cor_obj$estimate,3)
            pv <- signif(cor_obj$p.value,3)
            
            ggplot(ccdf, aes(xValue,yValue)) +
                geom_point(aes(colour=zValue))+
                xlab("Gene")+
                ylab("Cell")+
                ggtitle(paste0("Spearman: r = ",rv,", p value = ",pv))+
                labs(colour = "Cluster")+
                theme_bw()+ 
                theme(            
                    plot.title = element_text(hjust = 0.5),
                    axis.text = element_text(size = 12),
                    panel.grid = element_blank()
                )

        }    
    })  
    
    output$plot5 <- renderPlot({
        st <- input$st
        sampleName <- input$sampleName
        sc <- input$sc
        group <- input$group
        cellType <- input$cellType
        km <- as.character(input$km)

        if(is.list(cellMTX()) & is.list(clusterMTX()) & cellType%in%rownames(cellMTX()$visual.matrix))
        {
            olp <- intersect(colnames(cellMTX()$visual.matrix),colnames(clusterMTX()$visual.matrix))
            
            xValue <- clusterMTX()$visual.matrix[km,olp]
            yValue <- cellMTX()$visual.matrix[cellType,olp]
            
            ccdf <- data.frame(xValue=xValue,yValue=yValue)
            
            ggplot(ccdf, aes(xValue,yValue)) +
                geom_boxplot(aes(colour=xValue))+
                xlab("Cluster")+
                ylab("Cell")+
                labs(colour = "Cluster")+
                theme_bw()+ 
                theme(
                    axis.text = element_text(size = 12),
                    panel.grid = element_blank()
                )
        }    
    })  
    
    output$plot6 <- renderPlot({
        st <- input$st
        sampleName <- input$sampleName
        km <- as.character(input$km)

        if(is.list(averMTX()) & is.list(clusterMTX()) & "Average"%in%rownames(averMTX()$visual.matrix))
        {
            olp <- intersect(colnames(averMTX()$visual.matrix),colnames(clusterMTX()$visual.matrix))
            
            xValue <- clusterMTX()$visual.matrix[km,olp]
            yValue <- averMTX()$visual.matrix["Average",olp]
            
            ccdf <- data.frame(xValue=xValue,yValue=yValue)
            
            ggplot(ccdf, aes(xValue,yValue)) +
                geom_boxplot(aes(colour=xValue))+
                xlab("Cluster")+
                ylab("Gene")+
                labs(colour = "Cluster")+
                theme_bw()+ 
                theme(
                    axis.text = element_text(size = 12),
                    panel.grid = element_blank()
                )
        }    
    })  
    
})
