library(ggpubr)
library(viridis)

smokeStack <- function(df = cells, 
                       normBy = "zscore",
                       setX = F,
                       setY = F,
                       bins = 20,
                       binSize = F,
                       lims = F,
                       fixedCoord = T,
                       saveFile = "figures/coloredBar.pdf",
                       xdim = 5,
                       ydim = 0.5,
                       viridisPallet = "inferno"){
  if (setX == F){
    print(names(df))
    xName <- readline(prompt = "Which variable should be mapped on the x-axis? ")
  } else {
    xName <- setX
  }
  
  if (setY == F){
    varNumb <- as.numeric(readline(prompt = paste0("How many variables should be mapped across ", xName, "? ")))
    varMat <- data.frame("varNumb" = c(1:varNumb), "varName" = "HOLD")
    
    for (i in 1:varNumb){
      varMat[i,2] <- readline(prompt = paste0("What is variable number ", i, "? "))
    }
  } else {
    varMat <- data.frame("varNumb" = c(1:length(setY)), "varName" = setY)
  }
  
  if(binSize != F){
    if (is.numeric(binSize)){
      print(paste("Using bin size at", binSize))
    } else{
      print("Non-numeric bin size detected. Ignoring bin size")
      binSize <- F
    }
  }
  
  for (j in unique(varMat$varName)){
    if(xName %in% names(df) & j %in% names(df)){
      interim <- data.frame("Name" = df[,xName], "var" = df[,j])
      if (normBy == "zscore"){
        interim$var <- (interim$var-mean(interim$var))/sd(interim$var)
        legendName <- "Z-Score"
      } else if (normBy == "lq"){
        interim$var <- interim$var/quantile(interim$var)[2]
        legendName <- "% of LQ"
      } else if (normBy == "mean"){
        interim$var <- 100*interim$var/mean(interim$var)
        legendName <- "% of Mean"
      } else {
        legendName <- "Intensity"
      }
      if(is.double(lims())){
        interim <- interim[interim$Name >= lims[1] & interim$Name <= lims[2],]
      }
      interim$bin <- bins
      xNameBinRange <- (max(df[xName])-min(df[xName]))/bins
      binSet <- c(min(df[xName]))
      for (i in 1:bins-1){
        binSet <- append(binSet, rev(binSet)[1]+xNameBinRange)
      }
      binID <- 0
      for (i in 2:length(binSet)){
        binID <- binID+1
        if(nrow(interim[interim$Name >= binSet[i-1] & interim$Name <= binSet[i],])>0){
          interim[interim$Name >= binSet[i-1] & interim$Name <= binSet[i],]$bin <- binID
        }
      }
      
      interim$target <- j
      if(!exists("allData")){
        allData <- interim
      } else {
        allData <- rbind(allData, interim)
      }
    } else {
      print(paste0("Variable '", j, "' not detected in dataset. Please try again"))
    }
  }
  
  if (binSize != F & is.numeric(binSize)){
    for (i in unique(allData$bin)){
      ricky <- allData[allData$bin == i,]
      if (nrow(ricky) > binSize){
        ricky <- ricky[sample(nrow(ricky), binSize),]
      }
      if(!exists("ticky")){
        ticky <- ricky
      } else{
        ticky <- rbind(ticky, ricky)
      }
    }
    allData <- ticky
  }
  
  cBin <- c(1:bins)
  cBin <- cBin[c(T, F)]

  binSet <- binSet[2:length(binSet)]
  binSet <- round(binSet - (xNameBinRange/2),1)
  binSet <- binSet[c(T, F)]

  draft <- ggplot(data = allData, aes(x = bin, y = target, fill = var))+
    geom_tile()+
    theme_classic2()+
    scale_fill_viridis(option = viridisPallet)+
    xlab(xName)+
    ylab("")+
    guides(fill=guide_legend(title = legendName))+
    theme(legend.position = "bottom")+
    scale_x_continuous(breaks = cBin, labels = binSet)
  if(fixedCoord==T){
    draft  <- draft+coord_fixed()
  }
  print(draft)
  ggsave(saveFile)
}