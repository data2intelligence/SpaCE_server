library(ggplot2)

rotateDfr <- function(dfr, degree) {
  M <- dfr
  alpha <- pi * degree / 180 
  
  #rotation matrix
  rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2)
  
  #shift, rotate, shift back
  #M2 <- t(rotm %*% (t(M)-c(M[1,1],M[1,2]))+c(M[1,1],M[1,2]))
  M2 <- t(rotm %*% t(M))
  
  dfr <- data.frame(x=M2[,1],y=M2[,2])
  
  return(dfr)
}

normalizeInput <- function(item)
{
  item <- gsub(" +", " ", item, perl=T)
  item <- gsub("^ ", "", item, perl=T)
  item <- gsub(" $", "", item, perl=T)
  item <- unique(unlist(strsplit(item," ",perl=T)))
  toupper(item)
}

source("ST_data_new_path.R")
sts <- readLines("ST_data/name/dataSetName.txt")
scs <- readLines("ST_data/cellType/dataSetName.txt")
sts.y.reverse <- c("Stahl_2016_Breast.Cancer","Berglund_2018_Prostate.Cancer","Thrane_2018_Melanoma")

namePath <- "ST_data/name/"
imagePath <- "ST_data/image/"
preprocessDataLogCPMPath <- "ST_data/preprocessDataLogCPM/"

pointSize <- c(3,3,3,0.5,3)*2
names(pointSize) <- sts

figWidth <- 520
figHeight <- 470
