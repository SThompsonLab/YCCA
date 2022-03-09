#These functions are used to Core Analysis of Plaque IFA (CAP-IFA). After gating the dataset as infected or not infected (`variable`)
# you may all `foci_caller` to generate plaque IDs, a list of infected nearest neighbors within the "wiggle"-room, and the length of that list
# An option third function can be used to collect plaque data independent of the individual cell data

suppressPackageStartupMessages(library(hash))

#This is the first function which will generate the list of nearby infected cells for every cell within the dataset
# It is **IMPERATIVE** that the df variable is a dataset of cells from a single image. If your dataset has cells from multiple images, subset the dataset my image id and run the function on each subset separately
# The second function, which condenses the lists based on overlap, is automatically called
foci_caller <- function(df = cells,
                        variable = "nProtein",
                        variable_subset = "Positive",
                        midpointer_x = "XM_NUC_dna",
                        midpointer_y = "YM_NUC_dna", 
                        wiggle = 4,
                        fileName = "NN_all_cells.csv"){
  cat("Running foci_caller...\n")
  
  # Collect only the infected cells
  foci_cells <- df[df[variable] == variable_subset,]
  
  # Add some place-holder names
  df["focus"] <- 0
  df["NNs"] <- 0
  df["F_ids"] <- "NA"
  
  # Iterate through each infected cells to generate a list of distances from that cell to every cells in the image
  for (a in 1:nrow(foci_cells)){
    # Get the cell ID
    cellID_a <- foci_cells[a,]$Number
    # create a dummy dataset to work the math on
    distList <- df
    # Get the X and Y coordinates of the infected cell
    XM <- as.numeric(foci_cells[a,][midpointer_x])
    YM <- as.numeric(foci_cells[a,][midpointer_y])
    # Calculate the distances between the infected cell and every other cell and store it in the dummy dataset
    distList["distance"] <- sqrt((distList[midpointer_x]-XM)**2+(distList[midpointer_y]-YM)**2)
    # Order the dummy dataset  by the distances from closest to furthest
    distList <- distList[order(distList$distance),]
    # Generate a misses variable to use the wiggle variable
    misses <- 0
    focus <- distList[distList$Number == cellID_a,]

    # Now we go through the sorted dummy list to find a list of the closest, infected cells
    for (b in 1:nrow(distList[distList$Number != cellID_a,])){
      # Starting with itself, if the cell is infected, its number ID is added to the list
      if (distList[b,][variable] == variable_subset){
        # If the list doesn't exist, it is made and the number ID is added. Otherwise, its just added to the list
        focus <- rbind(focus, distList[b,])
        
        # If the next closest cell *isn't* infected, the misses variable is increased until it surpasses the wiggle variable
        # If this happens, the iteration through the sorted dummy dataset is ended and
        # the length and composition of the infected list is added to the original dataset
      } else{
        misses <- misses+1
        if (misses > wiggle){
          focus_size <- nrow(focus)
          df[df["Number"] == cellID_a,]$NNs <- focus_size
          
          F_id <- ""
          for (c in unique(focus$Number)){
            F_id <- paste(F_id, c, sep = " ")
          }
          F_id <- paste0(F_id, " ")
          df[df["Number"] == cellID_a,]$F_ids <- F_id
          rm(focus)
          break
        }
      }
    }
  }
  cat("foci called\n")
  # The new dataset is saved
  write.csv(df, fileName, row.names = F)
  # The next function will condense the lists of infected neighbors based on overlapping number IDs
  foci_compact(fileName = fileName,
               variable = variable,
               variable_subset = variable_subset)
}

foci_compact <- function(fileName = "assigned.csv", 
                         variable = "nProtein",
                         variable_subset = "Positive", 
                         overlap = 1){
  cat(paste0("Running foci_compact from ", fileName, ".\n"))
  cells_new <- read.csv(fileName)
  #First, we get the cells that are infected
  fList <- cells_new[cells_new[variable] == variable_subset,]
  #Then we'll order them from great to lowest
  fList <- fList[order(fList$NNs, decreasing = T),]
  # Then we create the hash to store the focal numbers and initialize the count
  focal_hash <- hash()
  foci_count <- 1
  cat("Hash generated...\n")
  
  for (a in 1:nrow(fList)){
    found <- F
    #Starting with the cell with the highest NNs value and working down, the NNs are split by the " " delimiter
    fIDs <- strsplit(as.character(fList[order(fList$NNs, decreasing = T),][a,]$F_ids), split = " ")
    fIDs <- fIDs[1][[1]][2:length(fIDs[1][[1]])]
    # If no hash entries are detected, the first infected cell NN becomes the first focus entry
    if (length(focal_hash) < 1){
      focal_hash[[as.character(foci_count)]] <- fIDs
      foci_count <- foci_count+1
      found <- T
      # Otherwise, the entries are compared the the previously entered hashes
    } else {
      hit <- 0
      for (b in fIDs){
        for (c in 1:length(focal_hash)){
          if(b %in% focal_hash[[as.character(c)]]){
            hit <- hit+1
            if(hit >= overlap){
              focal_hash[[as.character(c)]] <- unique(append(focal_hash[[as.character(c)]], fIDs))
              found <- T
              break
            }
          }
        }
      }
    }
    if (found == F){
      focal_hash[[as.character(foci_count)]] <- fIDs
      foci_count <- foci_count+1
    }
  }
  
  cat("Finding overlapping hashes...\n")
  for (a in ls(focal_hash)){
    found <- F
    b <- focal_hash[[a]]
    for (c in ls(focal_hash)){
      if (a != c){
        if (!is.na(as.numeric(table(b %in% focal_hash[[c]])["TRUE"]) >= overlap) & found == F){
          if (as.numeric(table(b %in% focal_hash[[c]])["TRUE"]) >= overlap){
            print(paste("matching", a, "with",c, "because they overlap by", as.numeric(table(b %in% focal_hash[[c]])["TRUE"])))
            focal_hash[[c]] <- unique(append(focal_hash[[c]], b))
            print(length(focal_hash[[c]]))
            del(a, focal_hash)
            found <- T
          }
        }
      }
      if (found==T){
        break
      }
    }
  }

  cat("Determining probabilities...\n")
  cells_new$focus <- 0
  cells_new$focal_size <- 0
  for (a in ls(focal_hash)){
    c <- length(unique(focal_hash[[a]]))
    for (b in focal_hash[[a]]){
      cells_new[cells_new$Number == b,]$focus <- a 
    }
    cells_new[cells_new["focus"] == a,]["focal_size"] <- c
  }
  cells_new$core_prob <- cells_new$NNs/cells_new$focal_size
  cells_new$core_prob_norm <- 0
  cells_new[is.na(cells_new$core_prob),]$core_prob <- 0
  for (a in unique(cells_new$focus)){
    cells_new[cells_new$focus == a,]$core_prob_norm <- 100*cells_new[cells_new$focus == a,]$core_prob/sum(cells_new[cells_new$focus == a,]$core_prob)
  }
  write.csv(cells_new, paste0("data/capIFA_",strsplit(fileName, "/")[[1]][2]), row.names = F)
  cat("Foci compacted.\n")
}

# This function should be run from the `files` directory and collects basic data on the plaques called in the earlier functions:
# Image ID, Infected status, plaque number in that image, number of cells in that plaque, the lowest and highest core probabilities, and the average X/Y coordinate of that plaque
plaque_data_collect <- function(){
  fileList <- list.files(path = "data/", pattern = "capIFA")
  if(exists("plaque_info")){
    rm(plaque_info)
  }
  
  if(exists("all_cells")){
    rm(all_cells)
  }
  
  for (a in fileList){
    print(paste("Opening file:",a))
    cells <- read.csv(paste0("data/", a))
    cells$top_core <- "Negative"
    for (b in unique(cells$focus)){
      cells[cells$focus == b & cells$core_prob ==  max(cells[cells$focus == b,]$core_prob),]$top_core <- "Positive"
    }
    if (!exists("all_cells")){
      all_cells <- cells
    } else {
      if(F %in% (names(cells) == names(all_cells))){
        for(i in names(all_cells)[!names(all_cells) %in% names(cells)]){
          cells[i] <- 0
        }
        for(i in names(cells)[!names(cells) %in% names(all_cells)]){
          all_cells[i] <- 0
        }
      }
      all_cells <- rbind(all_cells, cells)
    }
    for (c in unique(cells$focus)){
      tittering <- subset(cells, focus == c)
      if (!exists("plaque_info")){
        plaque_info <- data.frame("image" = unique(tittering$image),
                                  "time" = unique(tittering$time),
                                  "virus" = unique(tittering$virus),
                                  "infection" = unique(tittering$infection),
                                  "infected" = !unique(tittering$focus)==0,
                                  #"thirdStain" = unique(tittering$thirdStain),
                                  "plaque_number" = c,
                                  "plaque_size" = max(tittering$focal_size),
                                  "percent_G1" = 100*nrow(subset(tittering, ploidy == "2N" & edu == "Negative"))/nrow(tittering),
                                  "percent_G2" = 100*nrow(subset(tittering, ploidy == "4N" & edu == "Negative"))/nrow(tittering),
                                  "percent_S" = 100*nrow(subset(tittering, edu != "Negative"))/nrow(tittering),
                                  "low_prob" = min(tittering$core_prob),
                                  "high_prob" = max(tittering$core_prob),
                                  "Mean_x_position" = mean(tittering$X_NUC_dna),
                                  "Mean_y_position" = mean(tittering$Y_NUC_dna))
      } else {
        interim <- data.frame("image" = unique(tittering$image),
                              "time" = unique(tittering$time),
                              "virus" = unique(tittering$virus),
                              "infection" = unique(tittering$infection),
                              "infected" = !unique(tittering$focus)==0,
                              #"thirdStain" = unique(tittering$thirdStain),
                              "plaque_number" = c,
                              "plaque_size" = max(tittering$focal_size),
                              "percent_G1" = 100*nrow(subset(tittering, ploidy == "2N" & edu == "Negative"))/nrow(tittering),
                              "percent_G2" = 100*nrow(subset(tittering, ploidy == "4N" & edu == "Negative"))/nrow(tittering),
                              "percent_S" = 100*nrow(subset(tittering, edu != "Negative"))/nrow(tittering),
                              "low_prob" = min(tittering$core_prob),
                              "high_prob" = max(tittering$core_prob),
                              "Mean_x_position" = mean(tittering$X_NUC_dna),
                              "Mean_y_position" = mean(tittering$Y_NUC_dna))
        plaque_info <- rbind(plaque_info, interim)
      }
    }
  }
  write.csv(plaque_info, "data/plaque_data.csv", row.names = F)
  write.csv(all_cells, "data/experiment_data_w_coreProb_all_cells.csv", row.names = F)
}

