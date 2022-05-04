library(tidyverse)
library(ggpubr)

bud <- function(shape_file = "cona",
                max_distance = 20) {
  schema <- read_csv("schema.csv")
  setwd("files/")
  file_list <- list.files()
  file_list <- file_list[!grepl(".csv", file_list)]
  for (image_set in file_list){
    setwd(image_set)
    image_list <- list.files()
    image_list <- image_list[!grepl(".csv", image_list)]
    for (j in image_list){
      setwd(j)
      setwd("csvs/")
      cells <- read_csv(paste0(shape_file, ".csv"))
      names(cells)[1] <- "cell_id"
      
      cells$pRad <- cells$Perim./(2*pi)
      cells$aRad <- sqrt(cells$Area/pi)
      cells$radRat <- cells$pRad/cells$aRad
      cells$radRat <- cells$radRat/quantile(cells$radRat)[2]
      
      if ("dna.csv" %in% list.files()){
        dna <- read_csv("dna.csv")
        
        dna$pRad <- dna$Perim./(2*pi)
        dna$aRad <- sqrt(dna$Area/pi)
        dna$radRat <- dna$pRad/dna$aRad

        cells$dna_rois <- 0
        cells$total_dna_intDen <- 0
        cells$mean_dna_intDen <- 0
        cells$mean_dna_area <- 0

        for (i in 1:nrow(dna)){
          target_x <- dna[i,]$X
          target_y <- dna[i,]$Y
          target_intDen <- dna[i,]$IntDen
          target_area <- dna[i,]$Area
          
          interim <- cells
          interim$distance <- sqrt((cells$X-target_x)^2+(cells$Y-target_y)^2)
          minDist <- min(interim$distance)

          if(minDist < max_distance){
            cellID <- interim[interim$distance == minDist,]$cell_id[1]
            cells[cells$cell_id == cellID,]$dna_rois <- cells[cells$cell_id == cellID,]$dna_rois+1
            if (cells[cells$cell_id == cellID,]$mean_dna_intDen == 0){
              cells[cells$cell_id == cellID,]$mean_dna_intDen <- target_intDen
              cells[cells$cell_id == cellID,]$mean_dna_area <- target_area
            } else {
              cells[cells$cell_id == cellID,]$mean_dna_intDen <- (cells[cells$cell_id == cellID,]$mean_dna_intDen+target_intDen)/2
              cells[cells$cell_id == cellID,]$mean_dna_area <- (cells[cells$cell_id == cellID,]$mean_dna_area + target_area)/2
            }
            cells[cells$cell_id == cellID,]$total_dna_intDen <- cells[cells$cell_id == cellID,]$total_dna_intDen+target_intDen
          }
        }
      } else{
        print("DNA folder not detected, defaulting to size for cell cycle state...")
      }
      
      cells$log2_dna <- log(cells$mean_dna_intDen, 2)
      cells$image <- j
      cells$file <- image_set
      mData <- schema[schema$name_id == unique(cells$file),]
      for (mDatum in names(mData)[6:ncol(mData)]){
        cells[mDatum] <- mData[mDatum]
      }
      setwd("../")
      write.csv(cells, paste0(j, ".csv"), row.names = F)
      setwd("../")
    }
    setwd("../")
  }
  setwd("../")
}

#---------------------------------------------------

budCat <- function(unite = T){
  setwd("files/")
  dirList <- list.files()
  dirList <- dirList[!grepl(".csv", dirList)]
  for (dir in dirList){
    setwd(dir)
    imageList <- list.files()
    imageList <- imageList[!grepl(".csv", imageList)]
    for (image in imageList){
      setwd(image)
      if (!exists("totalList")){
        totalList <- read_csv(paste0(image, ".csv"))
      } else {
        interim <- read_csv(paste0(image, ".csv"))
        totalList <- rbind(totalList, interim)
      }
      setwd("../")
    }
    names(totalList)[1] <- "CellID"
    totalList$CellID <- paste0(totalList$image, "_", totalList$CellID)
    totalList <- totalList[,-2]
    write.csv(totalList, paste0(dir, ".csv"), row.names = F)
    if (unite == T){
      if(!exists("allCells")){
        allCells <- totalList
      } else {
        allCells <- rbind(allCells, totalList)
      }
    }
    rm(totalList)
    setwd("../")
  }
  setwd("../")
  if (unite == T){
    write.csv(allCells, "data/cells.csv", row.names = F)
  }
}

#------------------------------------------------------------

sell_cycles <- function(df = cells,
                        output_file = "data/cells_ycca.csv",
                        size_only = F,
                        re_graph = F,
                        re_graph_nulls = F,
                        G1_threshold = 1.15,
                        G2_threshold = 1.35,
                        max_distance = 100){
  draft <- ggplot(data = df, aes(x = radRat))+
    geom_density()
  print(draft+geom_vline(xintercept = G1_threshold, color = "red", linetype = 2, size = 2)+
          scale_x_continuous(breaks = c(0, 0.5, 1, 1.25, 1.5, 1.75, 2)))
  G1_check <- readline(prompt = paste0("Using ", G1_threshold, " as G1 threshold. Good? (Y/n): "))
  if (G1_check == "n"){
    G1_threshold <- as.numeric(readline(prompt = "What should the G1 threshold be: "))
    print(draft+geom_vline(xintercept = G1_threshold)+
            scale_x_continuous(breaks = c(0, 0.5, 1, 1.25, 1.5, 1.75, 2)))
    print(paste0("Using ", G1_threshold, " as G1 threshold..."))
  }
  df$state <- "G1"
  df[df$radRat > G1_threshold,]$state <- "not_G1"
  if (!"dna_rois" %in% names(df)){
    print("Not using DAPI files for cell cycle state...")
    print(draft+geom_vline(xintercept = G2_threshold, color = "red", linetype = 2, size = 2))
    G2_check <- readline(prompt = paste0("Using ", G2_threshold, " as G2 threshold. Good? (Y/n): "))
    if (G2_check == "n"){
      G2_threshold <- as.numeric(readline(prompt = "What should the threshold be: "))
      print(draft+geom_vline(xintercept = G2_threshold))
      print(paste0("Using ", G1_threshold, " as G2 threshold..."))
    }
    df[df$state != "G1" & df$radRat < G2_threshold,]$state <- "S"
    df[df$state != "G1" & df$radRat >= G2_threshold,]$state <- "G2"
  } else{
    df[df$state != "G1",]$state <- "S"
    df[df$state == "S" & df$dna_rois == 2,]$state <- "G2/M"
    if(nrow(df[df$dna_rois==0,]) > 0){
      df[df$dna_rois == 0,]$state <- "UNK"
    }
  }
  write.csv(df, output_file, row.names = F)
  
  if (re_graph){
    setwd("figures/")
    if (!"reMaps" %in% list.files()){
      dir.create("reMaps")
    }
    setwd("reMaps")
    for (fileType in unique(df$file)){
      if (!fileType %in% list.files()){
        dir.create(fileType)
      }
      setwd(fileType)
      for (image_number in unique(df$image)){
        if (re_graph_nulls){
          draft <- ggplot(data = df[df$image == image_number & df$file == fileType,], aes(x = XM, y = -YM, color = state))+
            geom_point(size = 2)+xlab("X position")+ylab("Y position")+theme_classic2()+
            theme(legend.position = "top")
        } else {
          draft <- ggplot(data = subset(df[df$image == image_number & df$file == fileType,], state != "UNK"), aes(x = XM, y = -YM, color = state))+
            geom_point(size = 2)+xlab("X position")+ylab("Y position")+theme_classic2()+
            theme(legend.position = "top")
        }
        print(draft)
        print(image_number)
        ggsave(paste0(image_number, ".png"), dpi = 300)
      }
      setwd("../")
    }
    setwd("../../")
  }
}

bud_explore <- function(df = cells,
                        ind_image = T,
                        output_file = "data/explored.csv"){
  scheme <- read_csv("schema.csv")
  johnny <- df
  johnny$dna_rois <- as.numeric(johnny$dna_rois)
  if (ind_image){
    johnny$target <- johnny$file
    johnny$file <- paste0(johnny$target, "_", johnny$image)
  }
  for (i in unique(johnny$file)){
    interim <- subset(johnny, file == i)
    if(!ind_image){
      mData <- scheme[scheme$name_id == i,]
    } else {
      mData <- scheme[scheme$name_id == unique(interim$target),]
    }
    
    if (!exists("expl")){
      expl <- data.frame("File" = i,
                         "Total" = nrow(interim),
                         "Mean_ConA_Area" = mean(interim$Area),
                         "Mean_ConA_Int" = mean(interim$Mean),
                         "Mean_ConA_IntDens" = mean(interim$RawIntDen),
                         "Mean_radRat" = mean(interim$radRat),
                         "Mean_DNA_rois" = mean(interim$dna_rois),
                         "Mean_DNA_Area" = mean(interim$mean_dna_area),
                         "Mean_DNA_Int" = mean(interim$mean_dna_intDen),
                         "Mean_DNA_IntDens" = mean(interim$total_dna_intDen),
                         "G0/G1" = 100*nrow(interim[interim$state == "G1",])/nrow(interim),
                         "S/G2" = 100*nrow(interim[interim$state == "S/G2",])/nrow(interim),
                         "G2/M" = 100*nrow(interim[interim$state == "G2/M",])/nrow(interim),
                         "Mitosis" = 100*nrow(interim[interim$state == "M",])/nrow(interim),
                         "UNK" = 100*nrow(interim[interim$state == "UNK",])/nrow(interim))
      for (mDatum in names(mData)[7:ncol(mData)-1]){
        expl[mDatum] <- mData[,mDatum]
      }
      if(ind_image){
        expl$target <- unique(interim$target)
      }
    } else {
      tick <- data.frame("File" = i,
                         "Total" = nrow(interim),
                         "Mean_ConA_Area" = mean(interim$Area),
                         "Mean_ConA_Int" = mean(interim$Mean),
                         "Mean_ConA_IntDens" = mean(interim$RawIntDen),
                         "Mean_radRat" = mean(interim$radRat),
                         "Mean_DNA_rois" = mean(interim$dna_rois),
                         "Mean_DNA_Area" = mean(interim$mean_dna_area),
                         "Mean_DNA_Int" = mean(interim$mean_dna_intDen),
                         "Mean_DNA_IntDens" = mean(interim$total_dna_intDen),
                         "G0/G1" = 100*nrow(interim[interim$state == "G1",])/nrow(interim),
                         "S/G2" = 100*nrow(interim[interim$state == "S/G2",])/nrow(interim),
                         "G2/M" = 100*nrow(interim[interim$state == "G2/M",])/nrow(interim),
                         "Mitosis" = 100*nrow(interim[interim$state == "M",])/nrow(interim),
                         "UNK" = 100*nrow(interim[interim$state == "UNK",])/nrow(interim))
      for (mDatum in names(mData)[7:ncol(mData)-1]){
        tick[mDatum] <- mData[,mDatum]
      }
      if(ind_image){
        tick$target <- unique(interim$target)
      }
      expl <- rbind(expl, tick)
    }
  }
  if (file.exists(output_file)){
    append_file <- readline(prompt = "Output file already exists. Append (Y/n)? ")
    if (append_file == "n"){
      output_file <- readline(prompt = "What should the new filename be: ")
      write_csv(expl, output_file)
    } else {
      interim <- read_csv(output_file)
      expl <- rbind(interim, expl)
      write_csv(expl, output_file)
    }
  } else {
    write_csv(expl, output_file)
  }
}

yeastCycleSelect <- function(df = cells,
                             default_radRat = 1.05,
                             default_log2dna = 10.5,
                             save_file = "data/cells_cells.csv",
                             explore = T,
                             explore_file = "data/cells_explored.csv"){
  df$deform <- F
  df$ploidy <- "2N"
  print("Gating on the radii deformation...")
  okay <- "n"
  iMag <- ggplot(data = df, aes(x = radRat))+
    geom_density(aes(color = dna_rois))+
    theme_classic2()
  print(iMag+geom_vline(color = "red", size = 1, linetype = 2, xintercept = default_radRat))
  print(paste0("Using a radii ratio threshold of: ", default_radRat))
  okay <- readline(prompt = "Does this look good (Y/n)? ")
  if (okay == "n"){
    while(okay == "n"){
      default_radRat <- as.numeric(readline(prompt = "What should the threshold be? "))
      print(paste0("Using a radii ratio threshold of: ", default_radRat))
      print(iMag+geom_vline(color = "red", size = 1, linetype = 2, xintercept = default_radRat))
      okay <- readline(prompt = "Does this look good (Y/n)? ")
    }
  } else {
    print(paste0("Using ", default_radRat, " as the threshold for radii deformations."))
  }
  df[df$radRat>default_radRat,]$deform <- T
  cat("\n")
  print("Gating on the DNA content now...")
  okay <- "n"
  iMag <- ggplot(data = df, aes(x = log2_dna))+
    geom_density(aes(color = dna_rois))+
    theme_classic2()
  print(iMag+geom_vline(color = "red", size = 1, linetype = 2, xintercept = default_log2dna))
  print(paste0("Using a ploidy threshold of: ", default_log2dna))
  okay <- readline(prompt = "Does this look good (Y/n)? ")
  if (okay == "n"){
    while(okay == "n"){
      default_log2dna <- as.numeric(readline(prompt = "What should the threshold be? "))
      print(paste0("Using a ploidy threshold of: ", default_log2dna))
      print(iMag+geom_vline(color = "red", size = 1, linetype = 2, xintercept = default_log2dna))
      okay <- readline(prompt = "Does this look good (Y/n)? ")
    }
  } else {
    print(paste0("Using ", default_log2dna, " as the threshold for DNA content."))
  }
  df[df$log2_dna > default_log2dna,]$ploidy <- "4N"
  df$state <- "G1"
  df[df$deform == F & df$ploidy == "4N",]$state <- "S/G2"
  df[df$deform == T & df$ploidy == "4N",]$state <- "G2/M"
  df[df$deform == T & df$ploidy == "2N",]$state <- "M"
  write.csv(df, save_file, row.names = F)
  if(explore==T){
    bud_explore(df, ind_image = F, output_file = explore_file)
  }
}

for (into in unique(cells$file)){
  intert <- subset(cells, file == into)
  yeastCycleSelect(df = intert, save_file = paste0("data/", into,".csv"))
}
