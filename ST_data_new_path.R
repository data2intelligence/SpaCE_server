sts <- readLines("ST_data/name/dataSetName.txt")

namePath <- "ST_data/name/"
rawDataPath <- "ST_data/rawData/"
preprocessDataPath <- "ST_data/preprocessData/"
preprocessDataLogCPMPath <- "ST_data/preprocessDataLogCPM/"
overviewPath <- "ST_data/overview/"
imagePath <- "ST_data/image/"
markerPath <- "ST_data/markerFig/"
clusteringPath <- "ST_data/clustering/"
deconvResPath <- "ST_data/deconvRes/"
deconvFigPath <- "ST_data/deconvFig/"

labSTPath <- "/data/Jiang_Lab/Data/ST/"
labSTPreprocessDataCPMPath <- "/data/Jiang_Lab/Data/ST/preprocessDataCPM/"

degrees <- list()
degrees[[sts[1]]] <- c(-26,-38,-170,160)
degrees[[sts[2]]] <- rep(0,12)
degrees[[sts[3]]] <- c(105,275,180,0,-66,25,-5,10)
degrees[[sts[4]]] <- rep(-90,2)
degrees[[sts[5]]] <- rep(0,2)

margins <- list()
for(i in 1:length(sts))
{
	sampleNames <- readLines(paste0(namePath,"sampleName_",sts[i],".txt"))
	names(degrees[[sts[i]]]) <- sampleNames
	
	if(i==1){
		margins[[sts[i]]][[sampleNames[1]]] <- c(0.2,0,0,0.5) #left, right, bottom, top
		margins[[sts[i]]][[sampleNames[2]]] <- c(0.2,0,0,0.5) 
		margins[[sts[i]]][[sampleNames[3]]] <- c(0.2,0,0,0.5) 
		margins[[sts[i]]][[sampleNames[4]]] <- c(0,0,0,0.5) 
	}else if (i==2){
		margins[[sts[i]]][[sampleNames[1]]] <- c(0.8,0.5,0.2,0.9)
		margins[[sts[i]]][[sampleNames[2]]] <- c(0.2,0.3,0.5,0.5)
		margins[[sts[i]]][[sampleNames[3]]] <- c(0,0.5,0.3,0.8)
		margins[[sts[i]]][[sampleNames[4]]] <- c(0,0,1,1)
		margins[[sts[i]]][[sampleNames[5]]] <- c(0,0.2,0.5,0.8)
		margins[[sts[i]]][[sampleNames[6]]] <- c(1.5,0.7,0.1,0.8)
		margins[[sts[i]]][[sampleNames[7]]] <- c(0.1,0.5,1,2)
		margins[[sts[i]]][[sampleNames[8]]] <- c(0.5,0,0.2,0.6)
		margins[[sts[i]]][[sampleNames[9]]] <- c(0.5,0,0.8,1.5)
		margins[[sts[i]]][[sampleNames[10]]] <- c(1.2,0.3,0.6,1)
		margins[[sts[i]]][[sampleNames[11]]] <- c(0.8,1,0,1)
		margins[[sts[i]]][[sampleNames[12]]] <- c(0,0,1,0.8)
	}else if (i==3){
		margins[[sts[i]]][[sampleNames[1]]] <- c(0,0,0.5,0.3) #left, right, bottom, top
		margins[[sts[i]]][[sampleNames[2]]] <- c(0.2,0,0.8,0.6) 
		margins[[sts[i]]][[sampleNames[3]]] <- c(0.3,1.3,1,1) 
		margins[[sts[i]]][[sampleNames[4]]] <- c(0.5,0.5,0.5,0.5) 
		margins[[sts[i]]][[sampleNames[5]]] <- c(0.5,0.5,0.5,0.6) #left, right, bottom, top
		margins[[sts[i]]][[sampleNames[6]]] <- c(0.2,0.2,1,1.3) 
		margins[[sts[i]]][[sampleNames[7]]] <- c(1,1,0.5,1.5) 
		margins[[sts[i]]][[sampleNames[8]]] <- c(1,0.8,1,1) 
	}else if (i==4){
		margins[[sts[i]]][[sampleNames[1]]] <- c(9,8,2.7,1.7) 
		margins[[sts[i]]][[sampleNames[2]]] <- c(12,11,4,1) 
	}else{
		margins[[sts[i]]][[sampleNames[1]]] <- c(0,1.8,0,0.5) 
		margins[[sts[i]]][[sampleNames[2]]] <- c(1,0,0,0)
	}
	
}

filterST <- function(st,sampleName,st.matrix.data)
{
	if(st=="Stahl_2016_Breast.Cancer")
	{
		if(sampleName=="Layer1_BC"){
			st.matrix.data <- st.matrix.data[,!colnames(st.matrix.data)%in%c("5.783x20.878","25.051x19.012","25.03x20.066","24.106x20.99","24.13x21.961","24.106x23.025","23.199x23.049","23.158x24.002","21.953x24.847","21.957x25.871")]
		}else if(sampleName=="Layer2_BC"){
			st.matrix.data <- st.matrix.data[,!colnames(st.matrix.data)%in%c("7.924x17.928","8.948x17.955","13.918x22.022","26.997x20.955","24.957x23.019","25.921x21.976","26.021x22.977","24.981x23.964","24.076x25.951")]
		}else if (sampleName=="Layer3_BC"){
			st.matrix.data <- st.matrix.data[,!colnames(st.matrix.data)%in%c("")]
		}else{
			st.matrix.data <- st.matrix.data[,!colnames(st.matrix.data)%in%c("12.083x9.046")]
		}
	}
	
	st.matrix.data
}




