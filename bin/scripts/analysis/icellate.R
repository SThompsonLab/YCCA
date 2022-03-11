suppressPackageStartupMessages(library(imager))
suppressPackageStartupMessages(library(magick))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

icellate <- function(targetCells,
                     folderName = "singleCells",
                     dMultiplier = 3,
                     verifySize = F,
                     fuse = F,
                     verifyImage = "overlay.png",
                     randomize = F,
                     samplingNumber = 5,
                     lineAnalyses = F,
                     numberOfLines = 2,
                     lineGraphs = F,
                     initialAngle = 0,
                     randomSamplerDenom = 10){
  #Set a parameter to allow size verification
  sizeVerified <- F
  fuseNames <- c()
  
  # First, we'll create the folder in case it doesn't exist, then go to it
  if(!"icellates" %in% list.files()){
    dir.create("icellates")
  }
  setwd("icellates")
  # Now, since each set of images are likely to represent *some* kind of cluster, we'll create of cluster folder thatn can be added to
  if(!folderName %in% list.files()){
    dir.create(folderName)
    
  }
  setwd(folderName)
  
  # If you want to randomly sample a certain number of cells from the population, this bit will do that
  if(randomize == T){
    targetCells <- targetCells[sample(nrow(targetCells), samplingNumber), ]
  }
  
  for (tCell in 1:nrow(targetCells)){
    targetCell <- targetCells[tCell,]
    #Now we get the target folder and the target image of the targetCell
    targetNameID <- unique(targetCell$name_id)
    targetImage <- unique(targetCell$image)
    
    # With this, we can make a special folder for JUST THIS CELL
    cellIdDir <- paste0(targetNameID, "-", unique(targetCell$CellID))
    if (!cellIdDir %in% list.files()){
      dir.create(cellIdDir)
      setwd(cellIdDir)
      # We can also now get a list of the images and move them into the new cellIdDir
      imageList <- list.files(path = paste0("../../../files/", targetNameID, "/", targetImage,"/PNGS"), 
                              pattern = ".png", full.names = T)
      for (i in imageList){
        file.copy(i, "./")
      }

      #Now to interate through the images, and them as a dataframe
      
      if (sizeVerified == F & verifySize == T){
        while (TRUE){
          interim <- load.image(verifyImage)
          interim <- as.data.frame(interim)
          # print("Image loaded")
          
          #First, lets figure out how big the area of the picture needs to be based on the imageType parameter
          areaName <- "Area"
          cellLength <- dMultiplier*round(sqrt(targetCell[areaName]/pi))
          cellLength <- as.integer(round(cellLength))
          
          # Now we identify the XY coordinate of the cell and trim the dataframe data to fit the calculated cellLength
          XName <- "X"
          XPos <- as.integer(round(targetCell[XName]))
          XMin <- round(XPos-cellLength)
          XMax <- round(XPos+cellLength)
          
          YName <- "Y"
          YPos <- as.integer(round(targetCell[YName]))
          YMin <- round(YPos-cellLength)
          YMax <- round(YPos+cellLength)
          
          # Now we remove all the pixels we don't care about
          interim <- interim[interim$x >= XMin & interim$x <= XMax,]
          interim <- interim[interim$y >= YMin & interim$y <= YMax,]
          
          #Now to reset the corrdinates
          interim$x <- 1+interim$x-min(interim$x)
          interim$y <- 1+interim$y-min(interim$y)
          
          # Now to convert the dataframe to an image object and save the image
          newPlot <- suppressWarnings(as.cimg(interim))
          plot(newPlot)
          print(paste0("dMultiplier is currently set to: ", dMultiplier))
          checkBox <- readline(prompt = "Does this look like the appropriate size (Y/n)? ")
          if (checkBox == "n"){
            dMultiplier <- as.numeric(readline(prompt = "What should the new dMultiplier be: "))
          } else{
            sizeVerified <- T
            break
          }
        }
      }
      
      # For each image...
      for (i in list.files()){
        # we load the image..
        interim_i <- load.image(i)
        # And convert it to a dataframe.
        interim <- as.data.frame(interim_i)
        # print("Image loaded")

        #First, lets figure out how big the area of the picture needs to be based on the imageType parameter
        areaName <- "Area"
        cellLength <- dMultiplier*round(sqrt(targetCell[areaName]/pi))
        cellLength <- as.integer(round(cellLength))
        # print("Cell length obtained")
        
        # Now we identify the XY coordinate of the cell and trim the dataframe data to fit the calculated cellLength
        XName <- "X"
        XPos <- as.integer(round(targetCell[XName]))
        XMin <- round(XPos-cellLength)
        XMax <- round(XPos+cellLength)
        # print("X position obtained")
        
        YName <- "Y"
        YPos <- as.integer(round(targetCell[YName]))
        YMin <- round(YPos-cellLength)
        YMax <- round(YPos+cellLength)
        # print("Y position obtained")

        # Now we remove all the pixels we don't care about
        interim <- interim[interim$x >= XMin & interim$x <= XMax,]
        interim <- interim[interim$y >= YMin & interim$y <= YMax,]
        # print("Pixels reduced")
        
        # Now we save the new centroid positions for possible line analyses
        newXPos <- 1+XPos-min(interim$x)
        newYPos <- 1+YPos-min(interim$y)

        #Now to reset the corrdinates
        interim$x <- 1+interim$x-min(interim$x)
        interim$y <- 1+interim$y-min(interim$y)
        
        # print("Centroids found")
        interim_now <<- interim
        
        # Now to convert the dataframe to an image object and save the image
        newPlot <- suppressWarnings(as.cimg(interim))
        # print("newPlot loaded")
        newPlot <- cimg2magick(newPlot, rotate = T)
        # print("newPlot rotated")
        image_write(newPlot, path = i, format = "png")
        # print("image saved")
        
        # If class averaging is called via the 'fuse' parameter AND the image is a square...
        if (fuse == T & max(interim$x) == max(interim$y)){
          # so we can call all the averaged images later, we'll make a list of each unique image
          if(!strsplit(i, ".png")[[1]][1] %in% fuseNames){
            fuseNames <- append(fuseNames, strsplit(i, ".png")[[1]][1])
          }
          # a new parameter is generated and assigned as the dataframe unless the variable already exists...
          if(!exists(paste0(strsplit(i, ".png")[[1]][1], "_fused"))){
            bidet <- interim
            bidet$number <- 1
            assign(paste0(strsplit(i, ".png")[[1]][1], "_fused"), bidet)
            assign(paste0(strsplit(i, ".png")[[1]][1], "_xOri"), newXPos)
            assign(paste0(strsplit(i, ".png")[[1]][1], "_yOri"), newYPos)
          } else {
            # in which case an interim dataframe is generated and the average is called in a declarable thing
            fuse_old <- get(paste0(strsplit(i, ".png")[[1]][1], "_fused"))
            # The number is increased
            fuse_old$number <- fuse_old$number+1
            fuse_interim <- interim
            
            #If the images are the same size, then the values are summed
            if(max(fuse_interim$x)==max(get(paste0(strsplit(i, ".png")[[1]][1], "_fused"))$x)){
              fuse_old$value <- fuse_old$value+fuse_interim$value
              assign(paste0(strsplit(i, ".png")[[1]][1], "_fused"), fuse_old)
              
            } else{
              # Now we use the midpoints XPos and YPos to align the dataframes

              if(nrow(fuse_old) > nrow(fuse_interim)){
                xShift <- get(paste0(strsplit(i, ".png")[[1]][1], "_xOri"))-newXPos
                yShift <- get(paste0(strsplit(i, ".png")[[1]][1], "_yOri"))-newYPos
                fuse_interim$x <- fuse_interim$x+xShift
                fuse_interim$y <- fuse_interim$y+yShift
              } else {
                xShift <- newXPos-get(paste0(strsplit(i, ".png")[[1]][1], "_xOri"))
                yShift <- newYPos-get(paste0(strsplit(i, ".png")[[1]][1], "_yOri"))
                fuse_old$x <- fuse_old$x+xShift
                fuse_old$y <- fuse_old$y+yShift
                assign(paste0(strsplit(i, ".png")[[1]][1], "_xOri"), newXPos)
                assign(paste0(strsplit(i, ".png")[[1]][1], "_yOri"), newYPos)
              }
              
              if(nrow(fuse_old[fuse_old$x == 0 | fuse_old$y == 0,])>=1){
                print(paste0("A total of", nrow(fuse_old[fuse_old$x == 0 | fuse_old$y == 0,]), "zero(s) was detected in the old"))
                fuse_old <- fuse_old[fuse_old$x > 0 & fuse_old$y > 0,]
              }
              if(nrow(fuse_interim[fuse_interim$x == 0 | fuse_interim$y == 0,])>=1){
                print(paste0("A total of ", nrow(fuse_interim[fuse_interim$x == 0 | fuse_interim$y == 0,]), "zero(s) was detected in the interim"))
                fuse_interim <- fuse_interim[fuse_interim$x > 0 & fuse_interim$y > 0,]
              }
              row.names(fuse_old) <- paste0(fuse_old$x, "_", fuse_old$y, "_", fuse_old$cc)
              row.names(fuse_interim) <- paste0(fuse_interim$x, "_", fuse_interim$y, "_", fuse_interim$cc)
              
              fuse_new <- merge(fuse_old[,1:4], fuse_interim, by = 0, all = T)
              fuse_new[is.na(fuse_new)] <- 0
              the_numbers <- as.numeric(unique(fuse_old$number))
              
              if(nrow(fuse_old) > nrow(fuse_interim)){
                fuse_export <- data.frame("x" = fuse_new$x.x,
                                          "y" = fuse_new$y.x,
                                          "cc" = fuse_new$cc.x,
                                          "value" = fuse_new$value.y+fuse_new$value.x,
                                          "number" = the_numbers)
              }else {
                fuse_export <- data.frame("x" = fuse_new$x.y,
                                          "y" = fuse_new$y.y,
                                          "cc" = fuse_new$cc.y,
                                          "value" = fuse_new$value.y+fuse_new$value.x,
                                          "number" = the_numbers)
              }
              #Now we create the new, averaged image
              assign(paste0(strsplit(i, ".png")[[1]][1], "_fused"), fuse_export)
            }
          }
        }
      }
      
      # Now to run the line analyses if they're supposed to be
      if(lineAnalyses == T){
        lineAnalysis(midX = newXPos,
                     midY = newYPos,
                     imageName = interim,
                     linesNumber = numberOfLines,
                     lineLength = cellLength,
                     graphing = lineGraphs,
                     offset = initialAngle, 
                     totalCorrPixFraction = randomSamplerDenom)
      }
      setwd("../")
    } else {
      print(paste0("Image folder exists, skipping ", cellIdDir))
    }
  }
  if(fuse==T){
    fuse_data <- data.frame("image" = fuseNames, "number" = 0)
    #this is where you call the assigned things and save them as images
    for (fused in fuseNames){
      # for each image name that was found, we get the image
      fused_image <- get(paste0(fused, "_fused"))
      #The intensities get rounded and then we drop the image numbers after saving that data
      fused_image$value <- fused_image$value/fused_image$number
      fuse_data[fuse_data$image == fused,]$number <- max(fused_image$number)
      fused_image <- fused_image[,1:4]
      # The averaged image is then saved as a png
      fusedPlot <- suppressWarnings(as.cimg(fused_image))
      fusedPlot <- cimg2magick(fusedPlot, rotate = T)
      image_write(fusedPlot, path = paste0(fused, ".png"), format = "png")
    }
    # The number of images used is then saved
    write.csv(fuse_data, "ClassAveragingMetrics.csv", row.names=F)
  }
  setwd("../../")
}

lineAnalysis <- function(midX,
                        midY,
                        imageName,
                        linesNumber = 2,
                        lineLength,
                        totalCorrPixFraction = 10, 
                        graphing = T, 
                        offset = 0){
  # First we make a directory to store the graphs and lines
  if (!"lineAnalyses" %in% list.files()){
    dir.create("lineAnalyses")
  }
  setwd("lineAnalyses")
  if (!"csv" %in% list.files()){
    dir.create("csv")
  }
  setwd("../")
  
  # Now we'll pull out the images images to combine them into a single pixel dataframe
  pngList <- list.files(pattern = ".png")
  for (pngImage in pngList){
    # We'll ignore the overlay image because we may be using more than RBG
    if (pngImage != "overlay.png"){
      # Load the image as an object
      lineImage <- load.image(pngImage)
      #Convert it to grayscale dataframe
      lineImage <- as.data.frame(grayscale(lineImage))
      lineImage$value <- lineImage$value-min(lineImage$value)
      # Rename the pixel value column as the stain ID
      names(lineImage)[3] <- strsplit(pngImage, ".png")[[1]][1]
      
      #Join the iterated extractions into a single dataset
      if (!exists("lineData")){
        lineData <- lineImage
      } else {
        lineData <- cbind(lineData, lineImage[strsplit(pngImage, ".png")[[1]][1]])
      }
      # Then we randomly select a fraction of the pixels to save for whole cell analysis
      randomPoints <- lineData[sample(nrow(lineData), nrow(lineData)/totalCorrPixFraction),]
      write.csv(randomPoints, "lineAnalyses/WholeCellRandom.csv", row.names = F)
    }
  }
      
  # now we pull out lines of pixel data that run through the center of the nucleus
      
  # If the linesNumber variable is 2, it just takes the x and y midpoint lines, 
  # otherwise it rotates radially though the center
  if (linesNumber < 3){
    maxX <- midX+(lineLength)
    minX <- midX-(lineLength)
    xLine <- subset(lineData, y == midY & x <= maxX & x >= minX)
    write.csv(xLine, file = "lineAnalyses/csv/xLine.csv", row.names = F)
    maxY <- midY+(lineLength)
    minY <- midY-(lineLength)
    yLine <- subset(lineData, x == midX & y <= maxY & y >= minY)
    write.csv(yLine, file = "lineAnalyses/csv/yLine.csv", row.names = F)
  } else {
    # we'll establish the initial angle as 0
    lineAngle <- 0
    # then we figure out the angle between lines
    angleNumber <- round(360/linesNumber)
    # now we iterate through each line to extract the line data
    for (line in 1:linesNumber){
      # We add the angle difference from the previous line and convert it radians
      lineAngle <- lineAngle+angleNumber
      angleRadians <- lineAngle*pi/180
      # Using the length parameter, we'll find the starting x and y pixels for the line
      xStart <- round(sin(angleRadians)*lineLength+midX)
      yStart <- round(cos(angleRadians)*lineLength+midY)
      xEnd <- round(midX-sin(angleRadians)*lineLength)
      yEnd <- round(midY-cos(angleRadians)*lineLength)
      
      # Since the change in x and y positions from the midpoint to the ending point mirrors the start we can determine the ending point
      
      # Now that we've got a mid-point and a starting point, we can get the slope of the line
      if (midX != xStart & midY != yStart){
        lineSlope <- -(midY-yStart)/(midX-xStart)
        yInt <- midY+(lineSlope*midX)
        newLine <- subset(lineData, 
                          x <= (max(xStart, xEnd)) &
                          x >= (min(xStart, xEnd)) &
                          y == -round(lineSlope*x-yInt))
      } else {
        newLine <- subset(lineData, 
                          x <= (max(xStart, xEnd)) &
                            x >= (min(xStart, xEnd)) &
                            y <= (max(yStart, yEnd)) &
                            y >= (min(yStart, yEnd)))
      }
      # Using the slope and the starting/ending points, we can trim the data to only pick the pixels along a line
      write.csv(newLine, file = paste0("lineAnalyses/csv/lineAngle_", lineAngle, ".csv"), row.names = F)
    }
  }
  
  # Now that we've got the line data, we can graph them because - why not?
  if (graphing == T){
    setwd("lineAnalyses/csv")
    lineFiles <- list.files(pattern = ".csv")
    for (lineFile in lineFiles){
      lineData <- read.csv(lineFile)
      if(!"pix" %in% names(lineData)){
        lineData$pix <- 1:nrow(lineData)
        colLength <- ncol(lineData)-1
        lineData <- lineData %>% pivot_longer(cols = 3:colLength, 
                                              names_to = "color",
                                              values_to = "intensity")
        lineData <- lineData[,c(3, 1, 2, 4, 5)]
        write.csv(lineData, file = paste0("longer_", lineFile), row.names = F)
      }
      linera <- ggplot(data = lineData, 
                       aes(x=pix, y = intensity, color=color))+
        geom_step(size = 1.5)+
        xlab("Relative pixel position")+
        ylab("Relative pixel intentsity")+
        theme_classic()+
        theme(axis.line = element_line(colour = "black", size = 2),
              axis.text = element_text(face = "bold", color = "black", size = 12),
              axis.title = element_text(face = "bold", color = "black", size = 20),
              legend.position = "top")
      print(linera)
      pdfName <- strsplit(lineFile, ".csv")[[1]][1]
      ggsave(filename = paste0("../", pdfName, ".pdf"))
    }
    setwd("../../")
  }
}