# #!/usr/bin/Rscript
# 
# #This assigns the arguments
# args = commandArgs(trailingOnly=TRUE)

# here's the actual function to automatically generate the directories based on the schema file
dirGen <- function(idType="name",
                   singleImage="auto",
                   imageNumber=5,
                   numberOfTargets=4){
  #First, lets generate the appropriate folders
  dirList <- list.files()
  if(!"files" %in% dirList){
    dir.create("files")
  }
  if(!"figures" %in% dirList){
    dir.create("figures")
  }
  if(!"data" %in% dirList){
    dir.create("data")
  }
  
  #reads the schema file
  schema <- read.csv("schema.csv")
  startingPoint <- numberOfTargets+2
  # if the passed argument is "number", the notebookID is used, otherwise a descriptive name is generated from the optional data tags
  if (idType == "number" | ncol(schema) < startingPoint){
    xList <- unique(schema$notebook_id)
  } else {
    for (a in startingPoint:ncol(schema)){
      if (exists("xList") & names(schema)[a] != "name_id"){
        xList <- paste(xList, schema[,a], sep = "_")
      } else if (names(schema)[a] != "name_id"){
        xList <- schema[,a]
      }
      schema$name_id <- xList
      write.csv(schema, "schema.csv", row.names = F)
    }
  }
  # Directories are touched, unless they already exist
  for(a in xList){
    dirName <- paste0("files/", a)
    if(!dir.exists(dirName))
    dir.create(dirName)
    if (singleImage == "single"){
      tip <- 1
      for (b in 1:imageNumber){
        imageDir <- paste0("image_", tip)
        dir.create(paste0("files/", a, "/", imageDir))
        tip <- tip+1
      }
    }
  }
}

runType <- readline(prompt = "Generating directories. Should this be done by name or number (NAME/number): ")
singleImages <- readline(prompt = "Are you taking single images, or automated capture (AUTO/single): ")
if (singleImages == "single"){
  imagesNumber <- readline(prompt = "How many images per condition (Default = 5): ")
  if (is.na(as.numeric(imagesNumber))) {
    imagesNumber <- 5
    print(paste0("Generating ", imagesNumber, " image folders in each condition..."))
  } else {
    imagesNumber <- as.numeric(imagesNumber)
    print(paste0("Generating ", imagesNumber, " image folders in each condition..."))
  }

}
ticko <- readline(prompt = "How many targets are listed in the schema (Default 4): ")
if (is.character(ticko)){
  ticko <- 4
} else {
  ticko <- as.numeric(ticko)
}
dirGen(idType = runType,
       singleImage = singleImages,
       imageNumber = imagesNumber,
       numberOfTargets = ticko)