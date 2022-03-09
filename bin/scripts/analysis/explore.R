explore <- function(fileName = "data/experiment.csv",
                    sortBy = "name_id",
                    categories = T,
                    cellCycle = T,
                    scheme = "schema.csv",
                    icellate = F,
                    vSize = F,
                    random = T,
                    extraMetrics = F,
                    saveFile = "data/experiment_explored.csv"){
  
  #First, filename is a dataset, it is used, otherwise it is opened as a csv. If the file is not a csv, an error is printed
  if(is.data.frame(fileName)){
    dFrame <- fileName
  } else if (grepl(".csv", fileName)){
    dFrame <- read.csv(file = fileName)
  } else {
    print("File or dataset not recognized.")
  }
  
  #A new, temporary copy of the name_id is created for seom reason?
  dFrame["placeHolder"] <- dFrame[sortBy]
  
  # if categories is true, you are asked home many categories to separate the data by
  if(categories == T){
    catNum <- as.numeric(readline(prompt = "How many categories would you like to explore: "))
    if (catNum > 0){
      # The names to choose from are listed
      print(names(dFrame)[1:length(dFrame)-1])
      # A temporary categories table is generated of an appropriate length
      catTab <- data.frame("catNum"=c(1:catNum))
      # An image_id name is generated for reasons?
      catTab$catName <- "imageID"

      # For each category number fiven, you're asked which variable should be used
      for (i in 1:catNum){
        catName <- readline(prompt = paste0("What is the name of category ", i, "? "))
        # After the category is generated, each row's placeholder variable gets update based on its value of the category
        dFrame$placeHolder <- paste0(dFrame$placeHolder, "_", catName, ".")
        dFrame$placeHolder <- paste0(dFrame$placeHolder, unlist(dFrame[catName]))
        catTab[catTab$catNum == i,]$catName <- catName
      }
      # Now that we have groups of unique category values, we cahnge the sortBy to use the lists rather than a name_id
      sortBy <- "placeHolder"
    }
  } else {
    catNum <- 0
  }
  
  # a new dataframe is generated to store all the frequency data by using the unique values of the placeholder variable
  datum <- data.frame("name_id" = unique(unlist(dFrame[sortBy])))
  
  # We read the schema to get metadata so we don't have to get it later
  schema <- read.csv(scheme)
  # store the meta names
  metaNames <- names(schema)[6:(length(names(schema))-1)]
  # Add the metadata to the final dataset row by row
  for (i in metaNames){
    datum[i] <- "Holder"
    for(j in 1:nrow(datum)){
      target <- unique(unlist(dFrame[dFrame$placeHolder == datum["name_id"][j,],][i]))
      datum[i][j,] <- target
    }
  }
  
  # We create another id that is sum of all the metadata so we can pull total data later
  if (categories == T & catNum > 0){
    for (noi in names(datum)[2:length(names(datum))]){
      if (!exists("check")){
        check <- unlist(datum[noi])
      } else {
        check <- paste(check, unlist(datum[noi]), sep = "_")
      }
    }
    datum$totalName <- check
  }
  
  #This pulls the category data per population and adds the population's unique variable for that category
  if (categories == T & catNum > 0){
    for (b in catTab$catName){
      datum[b] <- "Holder"
      for (c in 1:nrow(datum)){
        datum[b][c,] <- strsplit(datum["name_id"][c,], paste0(b, "."))[[1]][2]
        datum[b][c,] <- strsplit(datum[b][c,], "_")[[1]][1]
      }
    }
  }
  
  #Now we pull the total values for each
  datum$total <- 0
  for (totalSet in unique(datum$name_id)){
    datum[datum$name_id == totalSet,]$total <- nrow(dFrame[dFrame["placeHolder"]==totalSet,])
    
    if (icellate != F & is.numeric(icellate)){
      source("bin/scripts/analysis/icellate.R")
      if(icellate >= nrow(dFrame[dFrame["placeHolder"]==totalSet,])){
        icellate(targetCells = dFrame[dFrame["placeHolder"]==totalSet,], 
                 folderName = totalSet, 
                 verifySize = vSize, 
                 samplingNumber = nrow(dFrame[dFrame["placeHolder"]==totalSet,]), 
                 randomize = random)
      } else {
        icellate(targetCells = dFrame[dFrame["placeHolder"]==totalSet,], 
                 folderName = totalSet, 
                 verifySize = vSize, 
                 samplingNumber = icellate, 
                 randomize = random)
      }
    }
  }
  if (categories == T & catNum > 0){
    datum$subsetTotal <- datum$total
    datum$total <- 0
    for (totalSet in unique(datum$totalName)){
      totalCount <- sum(subset(datum, totalName == totalSet)$subsetTotal)
      datum[datum$totalName == totalSet,]$total <- totalCount
    }
    datum$subsetPercent <- round(100*datum$subsetTotal/datum$total, 2)
  }

  
  if (cellCycle == T){
    datum$MDI <- NA
    datum$sPhase <- NA
    datum$G1 <- NA
    datum$G2 <- NA
    datum$re_sPhase <- NA
    datum$G3 <- NA
    datum$S_MFI <- NA
    datum$S_MDI <- NA
    datum$rS_MFI <- NA
    datum$rS_MDI <- NA
    for (a in unique(datum$name_id)){
      set <- dFrame[dFrame["placeHolder"]==a,]
      datum[datum$name_id ==a,]$MDI <- mean(set$dna_norm)
      datum[datum$name_id ==a,]$sPhase <- 100*nrow(set[set$edu == "Positive",])/nrow(set)
      datum[datum$name_id ==a,]$G1 <- 100*nrow(set[set$edu == "Negative" & set$ploidy == "2N",])/nrow(set)
      datum[datum$name_id ==a,]$G2 <- 100*nrow(set[set$edu == "Negative" & set$ploidy == "4N",])/nrow(set)
      datum[datum$name_id ==a,]$re_sPhase <- 100*nrow(set[set$edu == "Positive" & set$ploidy == ">4N",])/nrow(set)
      datum[datum$name_id ==a,]$G3 <- 100*nrow(set[set$edu == "Negative" & set$ploidy == ">4N",])/nrow(set)
      datum[datum$name_id ==a,]$S_MFI <- mean(set[set$edu == "Positive" & set$ploidy != ">4N",]$edu_norm)
      datum[datum$name_id ==a,]$S_MDI <- mean(set[set$edu == "Positive" & set$ploidy != ">4N",]$dna_norm)
      datum[datum$name_id ==a,]$rS_MFI <- mean(set[set$edu == "Positive" & set$ploidy == ">4N",]$edu_norm)
      datum[datum$name_id ==a,]$rS_MDI <- mean(set[set$edu == "Positive" & set$ploidy == ">4N",]$dna_norm)
    }
  }
  
  if (extraMetrics != F){
    print(names(dFrame))
    metricsNumber <- as.numeric(readline(prompt = "How many extra variables would you like to add: "))
    metricsTab <- data.frame("metricNumber" = c(1:metricsNumber), "metricName" = "holder")
    for (i in 1:nrow(metricsTab)){
      metricsTab[i,]$metricName <- readline(prompt = paste0("Which metric should be added (", i, "): "))
    }
    for (metric in unique(metricsTab$metricName)){
      datum[paste0("Mean_", metric)] <- NA
      datum[paste0("SD_", metric)] <- NA
      datum[paste0("GMean_", metric)] <- NA
      datum[paste0("Median_", metric)] <- NA
      for (a in unique(datum$name_id)){
        set <- dFrame[dFrame["placeHolder"]==a,]
        datum[datum$name_id ==a,][paste0("Mean_", metric)] <- mean(as.numeric(unlist(set[metric])))
        datum[datum$name_id ==a,][paste0("SD_", metric)] <- sd(as.numeric(unlist(set[metric])))
        datum[datum$name_id ==a,][paste0("GMean_", metric)] <- exp(mean(log(as.numeric(unlist(set[metric])))))
        datum[datum$name_id ==a,][paste0("Median_", metric)] <- median(as.numeric(unlist(set[metric])))
      }
    }
  }
  
  
  if (file.exists(saveFile)){
    concatAway <- readline(prompt = "Save file with same name detected. Concatenate? (Y/n): ")
    if(concatAway != "n"){
      datum_total <- read.csv(saveFile)
      datum <- rbind(datum_total, datum)
    } else {
      saveFile <- readline(prompt = "What should the name of the new file be: ")
    }
  }
  write.csv(datum, file = saveFile, row.names = F)
}