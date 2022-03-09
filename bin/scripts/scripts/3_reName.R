#!/usr/bin/Rscript

reName <- function(){
  schema <- read.csv("schema.csv")
  nameTargets <- names(schema)[2:5]
  xList <- list.files(path = "files/")
  for (a in xList){
    if (!grepl(".ijm", x = a)){
      aList <- list.files(path = paste0("files/", a))
      for (aL in aList){
        if(a %in% unique(schema$notebook_id)){
          hit <- subset(schema, notebook_id == a)[1,]
        } else if (a %in% unique(schema$name_id)){
          hit <- subset(schema, name_id == a)[1,]
        }
        yList <- list.files(path = paste0("files/", a,"/", aL), pattern = ".tif")
        for (b in yList){
          if(T %in% sapply(nameTargets, grepl, strsplit(b, ".tif")[[1]][1])){
            b_sub <- nameTargets[sapply(nameTargets, grepl, strsplit(b, ".tif")[[1]][1])]
            x <- unique(hit[b_sub][,1])
            x <- paste0(x, ".tif")
            if (a != x){
              file.rename(paste0("files/", a,"/", aL, "/", b), paste0("files/", a, "/", aL, "/", x))
            }
          } else if(grepl("Overlay.tif", b)) {
            file.rename(paste0("files/", a,"/", aL, "/", b), paste0("files/", a, "/", aL, "/" ,"overlay.tif"))
          }
        }
        if("PNGS" %in% list.files(path = paste0("files/", a, "/", aL))){
          yList <- list.files(path = paste0("files/", a,"/", aL, "/PNGS"), pattern = ".png")
          for (b in yList){
            if(T %in% sapply(nameTargets, grepl, strsplit(b, "png")[[1]][1])){
              b_sub <- nameTargets[sapply(nameTargets, grepl, strsplit(b, ".tif")[[1]][1])]
              x <- unique(hit[b_sub][,1])
              x <- paste0(x, ".png")
              if (a != x){
                file.rename(paste0("files/", a,"/", aL, "/PNGS/", b), paste0("files/", a, "/", aL, "/PNGS/", x))
              }
            } else if(grepl("Overlay.png", b)) {
              file.rename(paste0("files/", a,"/", aL, "/PNGS/", b), paste0("files/", a, "/", aL, "/PNGS/overlay.png"))
            }
          }
        }
      }
    }
  }
}

# If you notice that the naming scheme isn't quite right, use this function to reset the original file architecture
# Then fix the schema file, and re-run reName
reset_names <- function(){
  xList <- list.files(path = "files/")
  setwd("files")
  for (a in xList){
    if (!grepl(".ijm", x = a)){
      setwd(a)
      aList <- list.files()
      for (aL in aList){
        setwd(aL)
        check <- list.files()
        if ("originals" %in% check){
          doubleCheck <- list.files(path = "originals/")
          if (length(doubleCheck) > 0){
            tifList <- list.files(pattern = "tif")
            for (tif in tifList){
              file.remove(tif)
            }
            dirList <- list.files()
            for (dir in dirList){
              if (dir != "originals"){
                unlink(dir, recursive = T)
              }
            }
            oriFiles <- list.files(path = "originals/")
            for (ori in oriFiles){
              file.copy(from = paste0("originals/", ori), to = paste0("./", ori))
            }
            unlink("originals", recursive = T)
          }
        }
        setwd("../")
      }
    }
    setwd("../")
  }
  setwd("../")
}