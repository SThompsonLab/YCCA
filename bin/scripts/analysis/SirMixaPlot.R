#SirMixaPlot version 8.2 - Danger! High Voltage
#----------------------------------------------------------------------------------------------------------

#This program is designed to take a csv (comma separated values) file and generate scatterplots based on columns.
# Changes from version 8.1:
#   -Changed the defualt EdU neg setting on normalizer()
#   -Added a work-around for concater() should rbinding fail due to missing columns
#   -imaGen() now also creates a concatanated csv of all the ROIs in the WholeCell folder with the _ROI_ tag
#   -Added runPCA() function to take a dataframe and generate a PCA biplot; still in beta testing

#Changes from SirMixaPlot version 8.0:
#   -Added some comments and typos
#   -Added the reMap() function for re-creating a facsimile of the original image

#Changes from SirMixaPlot version 7.0:
#   -Takes data in from the YggData.imj macro
#   -imaGen() is now much bigger and generates combined CSVs for nuclear vs whole cell comparisons

#Things still that need to be coded:
#   -Exploratory data analyses function called explore() made


#------------------------------------------------------------------------------------------------------------

# Functions for your pleasure

#--------------------------------------------------------------------------------------
#joiner() connects the sirmixaplot() with grapho() and can be used to graph new plots without running through sirmixaplot(). Default dataset is cells but any dataset can be used
joiner <- function(df = cells){
  #newAxis is reset for the changeAxis() function
  newAxis <<- FALSE
  
  #cats the dataset columns
  cat("\n")
  cat("Data categories detected:")
  cat("\n")
  print(sort(names(df)))
  cat("\n")
  
  #This sets the x and y variables
  px<<-readline(prompt = "px - What parameter is the x-axis: ")
  xname<<-readline(prompt = "What is the name of the x axis: ")
  cat("\n")
  py<<-readline(prompt = "py - What parameter is the y-axis: ")
  yname<<-readline(prompt = "What is the name of the y axis: ")
  cat("Working...")
  cat("\n")
  
  #calls the graphing function
  grapho(df)
}

#grapho() actually makes the scatter plots
grapho <- function(df = cells, 
                   X = px, 
                   Y = py, 
                   Xn = xname,
                   Yn = yname){
  
  #backup dataframe is generated
  interim <<- cells

  cat("Data frame generated")
  cat("\n")
  
  #Initial scatter plot is generated
  qq<<-ggplot(data=df, aes_string(X, Y))+
    theme_linedraw()+
    theme(plot.title = element_text(size=40, face="bold"))+
    theme(axis.title = element_text(size=30, face="bold"))+
    geom_point(size=2)+
    scale_x_continuous(name=Xn)+
    scale_y_continuous(name=Yn)+
    theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())
  
  if(!exists("newAxis")){
    newAxis <<- FALSE
  }
  if(newAxis == TRUE){
    changeAxis(store[1], store[2], store[3], store[4])
  }
  
  #Program finishes up
  cat("Plot generated")
  cat("\n")
  print(qq)

  cat("Call 'modthequad(df)' to apply quadrant analysis")
  cat("\n")
  cat("Call 'cull()' to remove certain populations from original dataset (DANGER, but less than before)")
  cat("\n")
  cat("Call 'eGod(df)' to apply ergodic analysis of S phase cells")
  cat("\n")
  cat("Call 'color_denisty()' to add density gradient to qq plot")
  cat("\n")
  cat("Call 'changeAxis(min X, max, X, min Y, max Y)' to change the axes of qq plots")
  cat("\n")
  cat("Call 'parser(df, df$variable)' to create a gated population")
  cat("\n")
}

#----------------------------------------------------------------
#sirmixaplot() loads the csv file (filo) and will call joiner() to start graphing if run is TRUE. This may lead to negPopCorrect if the appropriate names aren't present
#To re-run the negPopCorrect, simply call adjust as TRUE
sirmixaplot <- function(filo, 
                        run=TRUE, 
                        adjust=TRUE){
  #The script, m'lord
  
  #This is the beginning of the readline prompts for the program. User inputs desired output conditionals. Also, certain columns are added for safety...
  thingee <<- filo
  cat(paste("Opening file: ", thingee, sep=""))
  cat("\n")
  
  #This checks if the file has the '_cells', which signals it has been previously processed. Otherwise it extracts the appropriate rows and creates a cells file
  if (grepl("_cells.csv", thingee, fixed = TRUE)){
    cells <<- read.table(file=thingee,header=TRUE, fileEncoding = "latin1", sep = ",")
  } else {
    cells <<- read.table(file=thingee,header=TRUE, fileEncoding = "latin1", sep = ",")
    cat("Creating a new file, sir or madam.")
    cat("\n")
    thingee <<- paste0(substr(thingee, 1, nchar(thingee)-4), "_cells.csv")
    cat(paste(thingee, "created."))
    for (i in 1:ncol(cells)){
      if (grepl("Mean_", names(cells)[i])){
        cells[paste0("I", names(cells)[i])] <<- cells[,i]*cells$Area
      }
    }
    cells$log2_dna <<- log(cells$IntDen_NUC_dna, 2)
    for (i in 1:ncol(cells)){
      if (grepl("Mean_", names(cells)[i])){
        cells[paste0("L", names(cells)[i])] <<- log(cells[,i], 10)
      }
    }
  }
  
  #writes the file
  write.csv(cells, file = thingee, row.names = FALSE)
  cat("\n")
  
  #Checks that some size correctiong has occured and begins making the graph
  if(run==TRUE){
    if(TRUE %in% str_detect(names(cells), "ALIMean_.")){
      joiner(cells)
    } else if(adjust==TRUE){
      joiner(cells)
      gate()
    }
  }
}

#----------------------------------------------------------------------------------
#This function creates quadrants and calculates the percentage in each quadrant
#   Geom_line function can be altered to change the aesthetics of the quadrant lines

modthequad <- function(xi = 1.5,
                        yi = 1.5,
                        df = cells,
                        saveit = F, 
                        pause = T, 
                        graphit = T){
  if (pause == T){
    xi<-readline(prompt = "What is the x-intercept: ")
    yi<-readline(prompt = "What is the y-intercept: ")
  }
  quad <- c("Q1", "Q2", "Q3", "Q4")
  todos <- nrow(df)
  qt1 <<- nrow(subset(df, df[px] < xi & df[py] > yi))
  qt2 <<- nrow(subset(df, df[px] > xi & df[py] > yi))
  qt3 <<- nrow(subset(df, df[px] < xi & df[py] < yi))
  qt4 <<- nrow(subset(df, df[px] > xi & df[py] < yi))
  PercentTotal <- c(format(round(qt1/todos*100, 2), nsmall = 2), format(round(qt2/todos*100, 2), nsmall = 2), format(round(qt3/todos*100, 2), nsmall = 2), format(round(qt4/todos*100, 2), nsmall = 2))
  quads <- data.frame(quad, PercentTotal)
  xmin <- min(df[px])
  xmax <- max(df[px])
  ymin <- min(df[py])
  ymax <- max(df[py])
  qqmod <<- qq
  qqmod <<- qqmod+
    geom_hline(yintercept = as.numeric(yi), color = "red", linetype = "dashed", size = 2)+
    geom_vline(xintercept = as.numeric(xi), color = "red", linetype = "dashed", size = 2)+
    annotate("text", x = xmin, y = ymax, label = paste0("Q1: ", as.character(format(round((qt1/todos*100), 2), nsmall = 2)), "%"))+
    annotate("text", x = xmax, y = ymax, label = paste0("Q2: ", as.character(format(round((qt2/todos*100), 2), nsmall = 2)), "%"))+
    annotate("text", x = xmin, y = ymin, label = paste0("Q3: ", as.character(format(round((qt3/todos*100), 2), nsmall = 2)), "%"))+
    annotate("text", x = xmax, y = ymin, label = paste0("Q4: ", as.character(format(round((qt4/todos*100), 2), nsmall = 2)), "%"))
  print(quads)
  if(graphit==T){
    print(qqmod+geom_density_2d())
    if(saveit==T){
      ggsave(paste0(px,"_by_",py,"_mod.pdf"))
    }
    cat("\n")
    cat("call 'qqmod' to see plot with quadrants.")
    cat("\n")
  }
}

#---------------------------------------------------------------------------
#This function will remove certain points from the cells dataframe.
cull <- function(){
  print(sort(names(cells)))
  para <<- readline(prompt = "Which parameter should be culled: ")
  bORs <<- readline(prompt = "Data with values (b)igger, (s)maller, or (=) to the target: ")
  vale <<- readline(prompt = "What is the target value: ")
  remover<<-c()
  if (bORs=="b"){
    for (i in 1:nrow(cells)){
      checkr <<- cells[i, para]
      if (checkr > as.numeric(vale)){
        remover <<- c(remover, i)
      }
    }
    cells<<-cells[-remover,]
  } else if (bORs=="s"){
    for (i in 1:nrow(cells)){
      checkr <<- cells[i, para]
      if (checkr < as.numeric(vale)){
        remover <<- c(remover, i)
      }
    }
    cells<<-cells[-remover,]
  } else if (bORs=="="){
    for (i in 1:nrow(cells)){
      checkr <<- cells[i, para]
      if (checkr == vale){
        remover <<- c(remover, i)
      }
    }
    cells<<-cells[-remover,]
  } else{
    grapho(cells)
  }
  grapho(cells)
  cat("\n")
  cat("Previous 'cells' has been saved as 'interim'")
}

#----------------------------------------------------------------------------
#This function performs ergodic analysis on the graph, but can currently only do it along the x-axis, with a y-axis cut-off value.
#eGod will take your scatterplot and perform ergodic analysis on an x-axis

eGod <- function(df=cells){
  rate <<- as.integer(readline(prompt = "What is the rate value: "))
  bind <<- as.integer(readline(prompt = "How many bins: "))
  bn <<- as.integer(readline(prompt = "Remove number of final bins: "))
  hig <<- max(df[px])
  low <<- min(df[px])
  divine <<- (hig-low)/(bind-1)
  yt <<- as.integer(readline(prompt = "What is your S phase cutoff: "))
  binner <<- matrix(nrow = nrow(df), ncol = 1)
  
  sPhase <<- 0
  binner <<- c()
  binSize <<- c()
  
  for (i in 1:(bind-bn)){
    binSize[i] <<- low+(i*divine)
  }
  for (i in 1:nrow(df)){
    if (df[i,py]>yt){
      binner[i] <<- floor((df[i,px]-low)/divine)+1
      sPhase <<- sPhase+1
    }else{
      binner[i] <<- 0
    }
  }
  tot <- nrow(df)
  
  CpB <<- c()
  for (i in 1:(bind-bn)){
    ghetto <- sum(binner == i)
    CpB[i] <<- ghetto
  }
  
  a <- log(2)/rate
  Ft <- sPhase/tot
  bRate <<- c()
  for (i in 1:length(CpB)){
    bRate[i] <<- a*((2-Ft)/(CpB[i]/sPhase))
  }
  binNum <- c()
  for (i in 1:length(CpB)){
    binNum[i] <- i
  }
  justice <<- data.frame(binNum, binSize, CpB, bRate)
  hh <<- ggplot(justice, aes(binSize))+
    geom_bar(aes(weight=bRate))+
    xlab("Size of binning variable")+
    ylab("Rate")
  print(justice)
  print(hh)
}

#----------------------------------------------------------------
#This function creates new dataframes of gated populations.
gate <- function(df=cells,
                 gateName="negPop"){
  print(qq+geom_density_2d())
  cat(paste0("Gating out the ", gateName, " population:"))
  cat("\n")
  x1 <- readline(prompt = "Which value of pX should the gating begin: ")
  x2 <- readline(prompt = "Which value of pX should the gating end: ")
  y1 <- readline(prompt = "Which value of py should the gating begin: ")
  y2 <- readline(prompt = "Which value of py should the gating end: ")
  xSet <- c(x1, x2)
  ySet <- c(y1, y2)
  
  workaround <- subset(df, get(px) > as.numeric(xSet[1]) & get(px) < as.numeric(xSet[2]))
  workaround <<- subset(workaround, get(py) > as.numeric(ySet[1]) & get(py) < as.numeric(ySet[2]))
  assign(gateName, workaround, envir = .GlobalEnv)
  
  negPop_correct()
}

#This function uses a true negative population to adjust the background
negPop_correct <- function(){
  x <<- strsplit(px, "_")
  x <- paste0(x[[1]][length(x[[1]])-1], "_", x[[1]][length(x[[1]])])
  y <- strsplit(py, "_")
  y <- paste0(y[[1]][length(y[[1]])-1], "_", y[[1]][length(y[[1]])])
  px_background <<- mean(negPop[names(negPop)==paste0("Mean_", x)][,1])
  py_background <<- mean(negPop[names(negPop)==paste0("Mean_", y)][,1])
  if (x != "dna"){
    cells[paste0("AIMean_", x)] <<- cells[paste0("IMean_", x)]-(cells$Area*px_background)
    cells[paste0("AIMean_", x)] <<- cells[paste0("AIMean_", x)] + abs(min(cells[paste0("AIMean_", x)]))+1
    cells[paste0("ALIMean_", x)] <<- log(cells[paste0("AIMean_", x)], 10)
  }
  if (y != "dna"){
    cells[paste0("AIMean_", y)] <<- cells[paste0("IMean_", y)]-(cells$Area*py_background)
    cells[paste0("AIMean_", y)] <<- cells[paste0("AIMean_", y)] + abs(min(cells[paste0("AIMean_", y)]))+1
    cells[paste0("ALIMean_", y)] <<- log(cells[paste0("AIMean_", y)], 10)
  }
  write.csv(cells, file = thingee, row.names = FALSE)
  if (y == "NUC_edu" | x == "NUC_edu"){
    normalizer()
  }
}

#----------------------------------------------------------------
#This function will automatically pseudo-color density of a cell cycle profile and provides the basic functions for density pseudo-coloring
get_density <- function(px, 
                        py,
                        ...){
  theme_set(theme_bw(base_size = 16))
  dens <<- MASS::kde2d(px, py, ...)
  ix <<- findInterval(px, dens$x)
  iy <<- findInterval(py, dens$y)
  ii <<- cbind(ix, iy)
  return(dens$z[ii])
}

color_density<-function(){
  cells$density <<- get_density(cells[px][,1], cells[py][,1], n=100)
  qq <<- ggplot(data = cells, aes_string(px, py)) + geom_point(aes(color = density)) + scale_color_viridis()+xlab(xname)+ylab(yname)
  
  if(newAxis == TRUE){
    changeAxis(store[1], store[2], store[3], store[4])
  }
  print(qq)
}

#----------------------------------------------------------------
#This function changes the axis of qq files via changeAxis(X0, X1, Y0, Y1)
changeAxis <- function(w,
                       x,
                       y,
                       z){
  newAxis <<- TRUE
  store <<- c(w,x,y,z)
  qq<<-qq+coord_cartesian(xlim = c(w,x), ylim = c(y, z))
  print(qq+geom_density_2d())
}

#----------------------------------------------------------------
#This function asks you to define a G1 peak and re-assigns that value to  '1'
normalizer <- function(){
  
  #create a column holder (just in case...)
  cells$edu<<-"holder"
  cells$ploidy<<-"holder"
  cat("Determining the G1 population now.")
  cat("\n")
  
  #This part calls an EdU negative popuation and creates the normalized values
  bigBin <- which.max(density(cells$ALIMean_NUC_edu)$y)
  edu_neg <- density(cells$ALIMean_NUC_edu)$x[bigBin]+0.5
  edu_hist <- ggplot(cells, aes(ALIMean_NUC_edu))+geom_density()+geom_vline(xintercept = edu_neg)+xlab("Log EdU")
  print(edu_hist)
  g1_good <- readline(prompt = paste0("Is ", format(round(edu_neg, 2), nsmall = 2), " representative of the EdU cutoff? (Y/n) "))
  if (g1_good == "n"){
    edu_neg = as.numeric(readline(prompt = "What value should EdU be cutoff as negative? "))
    edu_hist <- ggplot(cells, aes(ALIMean_NUC_edu))+geom_density()+geom_vline(xintercept = edu_neg)+xlab("Log EdU")
    print(edu_hist)
  }
  cat(paste0("Using ", format(round(edu_neg, 2), nsmall = 2), " as the EdU cutoff."))
  cat("\n")
  for (i in 1:nrow(cells)){
    if (cells$ALIMean_NUC_edu[i] > edu_neg){
      cells$edu[i] <<- "Positive"
    } else {
      cells$edu[i] <<- "Negative"
    }
  }
  eduNegCells <-subset(cells, ALIMean_NUC_edu < edu_neg)
  cells$edu_norm <<- (cells$ALIMean_NUC_edu+1)-mean(eduNegCells$ALIMean_NUC_edu)
  
  #This part calls the 2N peak that is EdU-negative   
  bigBin <- which.max(density(eduNegCells$log2_dna)$y)
  diploid <- density(eduNegCells$log2_dna)$x[bigBin]
  dna_hist <<- ggplot(eduNegCells, aes(log2_dna))+geom_density()+geom_vline(xintercept = diploid)+xlab("DNA content")
  print(dna_hist)
  g1_good <- readline(prompt = paste0("Is ", format(round(diploid, 2), nsmall = 2), " representative of the 2N peak? (Y/n) "))
  if (g1_good == "n"){
    diploid = as.numeric(readline(prompt = "What is the value of the 2N peak? "))
    dna_hist <<- ggplot(eduNegCells, aes(log2_dna))+geom_density()+geom_vline(xintercept = diploid)+xlab("DNA content")
    print(dna_hist)
  }
  cat(paste0("Using ", format(round(diploid, 2), nsmall = 2), " as the 2N peak value."))
  cat("\n")
  
  
  #This part bins each point's ploidy        
  cells$dna_norm <<- (cells$log2_dna+1)-diploid
  
  for (i in 1:nrow(cells)){
    if (cells$dna_norm[i] < 1.6) {
      cells$ploidy[i] <<- "2N"
    } else if (cells$dna_norm[i] > 2.6) {
      cells$ploidy[i] <<- ">4N"
    } else{
      cells$ploidy[i] <<- "4N"
    }
  }
  
  ckF <- readline(prompt = paste0("Is the filename ", thingee, "? (Y/n) "))
  if (ckF != "n"){
    write.csv(cells, file = thingee, row.names = FALSE)
  } else{
    thingee <<-readline(prompt = 'What is the name of the file: ')
    write.csv(cells, file = thingee, row.names = FALSE)
  }
  px<<-"dna_norm"
  xname <<- "DNA content (Log 2)"
  py <<- "edu_norm"
  yname <<- "EdU content (Log 10)"
  changeAxis(0.5,3.5,0.5,2.5)
  grapho(cells)
}

#----------------------------------------------------------------
#this function just makes a histogram for you (yay!)
makehisto<-function(df=cells, va="log2_dna"){
  hh<<-ggplot(df, aes(x = df[va]))+geom_density()+
    theme_classic()+
    theme(axis.line = element_line(color = "black", size = 1.5), 
          axis.ticks = element_line(color = "black", size = 1.5), 
          axis.text = element_text(size = 24, family = "sans", color = "black"),
          axis.title = element_text(size = 36, family = "sans", color = "black"))+
    xlab("Variable")
  print(hh)
}

#-----------------------------------------------------------------
parser <- function(df=cells, va){
  catNam <<- readline(prompt = "What is the name of this variable: ")
  catNum <<- as.integer(readline(prompt = "How many categories: "))
  cells$holder <<- "NA"
  catz <<- c()
  for (i in 1:catNum){
    catz[i]<-readline(prompt = paste0("What is the name of category ", i, " of ", catNam, ": "))
  }
  cat("\n")
  for (i in catz){
    owl <<- as.numeric(readline(prompt = paste0("What is the low value of ", i, ": ")))
    hige <<- as.numeric(readline(prompt = paste0("What is the high value of ", i, ": ")))
    for (j in 1:nrow(cells)){
      if(va[j] >= owl & va[j] <= hige){
        cells$holder[j]<<-i
      }
    }
    cat("\n")
  }
  if (catNam %in% colnames(cells)){
    cells[match(catNam, names(cells))] <<- NULL
    colnames(cells)[match("holder", names(cells))] <<- catNam
  } else{
    colnames(cells)[match("holder", names(cells))] <<- catNam
  }
  write.csv(cells, file = thingee, row.names = FALSE)
}

#-----------------------------------------------------------------
# So this is the compensating function in case (god forbid) you have spectral overlap. Make sure to run this having graphed the offending colors against eachother
# IMPORTANT - The skewed color must be graphed on the y axis
# IMPORTANT - This should be obsolete if purgo.R has been run correctly
compensate <- function(df = cells){
  # First it generates a datframe to hold point. This is kept open to eventaully morph this into a shape fitting gate
  gridIron <- data.frame(X=c(1), Y=c(1))
  cat("Thank you for choosing coompensate(). To begin, lets assign a first point. Be sure the intended line follows the overlap line:")
  # Then you pick a point on the line of the overlapping population
  gridIron$X[1] <- as.numeric(readline(prompt = "What is the x-value of the first point? "))
  gridIron$Y[1] <- as.numeric(readline(prompt = "What is the y-value of the first point? "))
  cat("\n")
  print(qq+geom_point(data = gridIron, aes(x=gridIron$X, y=gridIron$Y, color = "red", size = 16)))
  cat(paste0("Great, I have added that point. Now lets move on to the second point"))
  gridIron <- rbind(gridIron, c(1, 1))
  # Then you pick a second point
  gridIron$X[2] <- as.numeric(readline(prompt = "What is the x-value of the second point? "))
  gridIron$Y[2] <- as.numeric(readline(prompt = "What is the y-value of the second point? "))
  print(qq+geom_point(data = gridIron, aes(x=gridIron$X, y=gridIron$Y, color = "red", size = 16))+
          geom_segment(aes(x = gridIron[1,1], y = gridIron[1,2], xend = gridIron[2,1], yend = gridIron[2,2], color = "red", size=16)))
  # Asks if these two points lie on the 
  #print(gridIron)
  yORn <<- readline(prompt = "Does this look right (y/n)? ")
  if (yORn == "n"){
    cat("Sorry, lets try that again.")
    cat("\n")
    compensate()
  } else {
    cat("Excellent, generating the linear equation now:")
    cat("\n")
    rise <- gridIron[2,2]-gridIron[1,2]
    run <- gridIron[2,1]-gridIron[1,1]
    slop <- rise/run
    yint <- gridIron[1,2]-(gridIron[1,1]*slop)
    print(qq+geom_point(data = gridIron, aes(x=gridIron$X, y=gridIron$Y, color = "red", size = 16))+
            geom_abline(slope = slop, intercept = yint, size = 2, color = "blue"))
    cat(paste0("Slope has been determined to be: ", slop))
    cat("\n")
    cat(paste0("Y-intercept has been determined to be: ", yint))
    yORn <- readline(prompt = "Does this look right (y/n)? ")
    if (yORn == "n"){
      cat("Sorry, lets try that again.")
      cat("\n")
      compensate()
    } else {
      cat("Applying compensation.")
      check <- df
      check[py] <- check[py]-(slop*check[px]+yint)
      cat("Regraphing:")
      grapho(check)
      yORn <- readline(prompt = "Does this look right (y/n)? ")
      if (yORn == "n"){
        cat("Sorry, lets try that again.")
        cat("\n")
        compensate()
      } else {
        cells <<- check
        cat(paste0("Cells saved with compensated ", py))
        yORn <- readline(prompt = "Would you like to recalculate the means (y/n)? ")
        if (yORn == "y"){
          for (i in 1:ncol(cells)){
            if (grepl(py, names(cells)[i])){
              cells[paste0("I", names(cells)[i])] <<- cells[,i]*cells$Area
            }
          }
          for (i in 1:ncol(cells)){
            if (grepl(py, names(cells)[i])){
              cells[paste0("L", names(cells)[i])] <<- log(cells[,i]+1, 10)
            }
          }
          cat("Means recalculated. Thank you for your patience and patronage.")
          write.csv(cells, thingee, row.names = FALSE)
        }
      }
    }
  }
}

#------------------------------------------------------

# This function take the nuclear locations and sizes and re-maps the original image
#   df is the dataframe, defaulted to 'cells'
reMap <- function(df = cells, 
                  #   highlight is a vector containing the Number(s) to be highlighted in remapped image
                  highlight = F, 
                  #   X and Y refer to the beginning of the X and Y position names in the dataframe, defaulted to the nuclear geometric position
                  Y = "X_NUC_", 
                  X = "Y_NUC_", 
                  #   S is the size of the ROI, defaulted to nucelar size
                  S = "Area_NUC_", 
                  #   Xn and Yn are the axis names and don't really need to be changed
                  Xn = "X position", 
                  Yn = "Y position",
                  #   colz is a vector of colors to use for the highlights variable. Be sure to have enough colors for the number of populations in the highlights variable
                  colz = cb_black) {
  #These set the variable names using regex to anchor tot he front of the variable name
  X <- names(df)[str_detect(names(df), paste0("^",X))]
  Y <- names(df)[str_detect(names(df), paste0("^",Y))]
  S <- names(df)[str_detect(names(df), paste0("^",S))]
  #This creates a ggplot object for graphing
  test <- ggplot(data = df, aes(y = -unlist(df[X]), x = unlist(df[Y]), size = unlist(df[S])))
  #If the highlight variable is a vector and is not F, its graphed here
  if (is.vector(highlight) & highlight != F) {
    #print(1)
    check <- df
    check$highlight <- check$Number %in% highlight
    qq <<- test+geom_point(aes(color = check$highlight))+
      ylab(Yn)+
      xlab(Xn)+
      scale_color_manual(values = colz)+
      theme_classic()+
      theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())+
      theme(plot.title = element_text(size=40, face="bold"))+
      theme(axis.title = element_text(size=30, face="bold"))
  } 
  # If highlight is a not a vector (such as variable called) and isn't F, then its graphed
  else if (length(highlight) != 1 & highlight != F){
    #print(2)
    qq <<- test+geom_point(aes(color = highlight))+
      ylab(Yn)+
      xlab(Xn)+
      theme_classic()+
      theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())+
      theme(plot.title = element_text(size=40, face="bold"))+
      theme(axis.title = element_text(size=30, face="bold"))+
      scale_color_manual(values = colz)
  }
  # Finally, if highlight IS FALSE, then its graphed without a color option
  else {
    #print(3)
    qq <<- test+geom_point()+
      ylab(Yn)+
      xlab(Xn)+
      theme_classic()+
      theme(axis.text = element_text(face='bold', size=16), axis.ticks = element_blank())+
      theme(plot.title = element_text(size=40, face="bold"))+
      theme(axis.title = element_text(size=30, face="bold"))
  }
  #The plot is graphed
  print(qq)
}

#-----------------------------------------------------------------
# Runs either a PCA of t-SNE
runCA <- function(df=cells,
                   plotType = "tsne",
                   saveFile = "figures/pcaPlot.png",
                   subz = F, 
                   groupz = T){  
  cellist <- df
  
  # Get rid of NAs and Infs
  cellist[is.na(cellist)] <- 0
  cellist[cellist == Inf] <- 0
  cellist[cellist == -Inf] <- 0 

  if (groupz ==  T){
    elly <- T
    print(names(cellist))
    grope <- readline(prompt = paste0("What should the ", plotType, " be grouped by from the 'cellist' dataframe: "))
    groupList <- unlist(df[grope])
  }
  # get rid of non-numerics
  num <- unlist(lapply(df, is.numeric))
  cellist_nums <- df[ , num]
  
  # Get only the nuclear columns & remove NA
  if (subz != F){
    cellist_nums <- cellist_nums[str_detect(names(cellist_nums), subz) | str_detect(names(cellist_nums), "Number")]
    cellist <- cellist[str_detect(names(cellist), subz) | str_detect(names(cellist), "Number")]
  }
  
  #Get rid of non-unique rows
  cellist <- unique(cellist)
  cellist_nums <- unique(cellist_nums)
  
  if (plotType != "pca" & plotType != "tsne"){
    plotType <- readline(prompt = "Are you trying to make a tsne or pca plot: ")
  }
  
  # Make the pca object
  if(plotType == "pca"){
    # Get rid of columns with no variance
    x <- foo(cellist_nums)
    cellist_nums <- cellist_nums[,-x]
    cellist_nums[is.na(cellist_nums)] <- 0
    
    cell.pca <- prcomp(cellist_nums, center = TRUE, scale = TRUE)
    percentage <- round(cell.pca$sdev / sum(cell.pca$sdev) * 100, 2)

    if (groupz == T){
      pca_out <- data.frame(x=cell.pca$x[,1], y = cell.pca$x[,2], groupBy = cellist[,grope])
      pp <<- ggplot(data = pca_out, aes(x = x, y = y))+
        geom_point(aes(color = groupBy), size = 3, alpha = 0.5)+
        geom_density_2d(aes(color = groupBy))+
        xlab(paste0("PC1 (", percentage[1], " %)"))+
        ylab(paste0("PC2 (", percentage[2], " %)"))+
        theme_classic2()
    } else {
      pca_out <- data.frame(x=cell.pca$x[,1], y = cell.pca$x[,2])
      pp <<- ggplot(data = pca_out, aes(x = x, y = y))+
        geom_point(size = 3, alpha = 0.5)+
        geom_density_2d()+
        xlab(paste0("PC1 (", percentage[1], " %)"))+
        ylab(paste0("PC2 (", percentage[2], " %)"))+
        theme_classic2()
    }
  } else{
    cellist_m <- as.matrix(cellist_nums)
    cellist_m[is.na(cellist_m)] <- 0
    tsne_out <- Rtsne::Rtsne(cellist_m, perplexity = 100)
    if (groupz == T){
      tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], groupBy = cellist[,grope])
      pp <<- ggplot(data = tsne_plot, aes(x=x, y=y))+geom_point(aes(color = groupBy), size = 3, alpha = 0.5)+geom_density_2d(aes(color = groupBy))+theme_classic2()
    } else{
      tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2])
      pp <<- ggplot(data = tsne_plot, aes(x=x, y=y))+geom_point(size = 3, alpha = 0.5)+geom_density_2d()+theme_classic2()
    }
  }
  print(pp)

  #Report how much of the data is shown
  howManyLeft <- round((nrow(cellist_nums)*ncol(cellist_nums))/(nrow(df)*ncol(df))*100, digits = 1)
  cat(paste0("Approximately ", howManyLeft, "% of the data is represented after trimming."))
  cat("\n")
  
  
  ggsave(saveFile, dpi = 300, height = 4, width = 4)
}

foo <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(!out > 1)
  unlist(want)
}

#-----------------------------------------------------------------
#nearest neighbor function for determining the nearest cell from one population to any cell of another population
nearestN <- function(home=cells, target=cells, label="CellToCell"){
  cells[paste0("nearest_",label)] <<- "NotHome"
  cells[paste0("nearest_",label,"_distance")] <<- 0
  for (i in 1:nrow(cells)){
    if(cells$Number[i] %in% unique(home$Number)){
      pool <- subset(target, Number != cells$Number[i])
      pool$distance <- sqrt((pool$X_NUC_dna-cells$X_NUC_dna[i])**2+(pool$Y_NUC_dna-cells$Y_NUC_dna[i])**2)
      cells[paste0("nearest_",label)][i,] <<- subset(pool, distance == min(pool$distance))$Number[1]
      cells[paste0("nearest_",label, "_distance")][i,] <<- subset(pool, distance == min(pool$distance))$distance[1]
    }
  }
}

#-----------------------------------------------------------------
# This function uses the general principle of nearestN but assumes you're clustering cells based on a gated variable.
# The function sorts all other cells based on their distance to each cell, then counts and records the nearby cells that are within a selected category until it encounters a nearby cell that is within the edge category
# It then stops can records the size of the cluster centered around that cell, and the numberIDs of those cells within the cluster
# After characterizing the clusters of each cell, it then compares the various clusters
clusterCount <- function(df = cells, clustPop="tag", center="Positive", edge="Negative", clusterCutoff=0.5){
  cells["clusterSize"] <<- 0
  cells["clusterMates"] <<- "None"
  cells["clusterID"] <<- "None"
  clusterID <- 1
  home <- cells[cells[clustPop]==center,]
  for (i in 1:nrow(cells)){
    if(cells$Number[i] %in% unique(home$Number)){
      pool <- subset(cells, Number != cells$Number[i])
      pool$distance <- sqrt((pool$X_NUC_dna-cells$X_NUC_dna[i])**2+(pool$Y_NUC_dna-cells$Y_NUC_dna[i])**2)
      pool <- pool[order(pool$distance),]
      clusterSize <- 1
      clusterMates <- c(cells$Number[i])
      for (j in 1:nrow(pool)){
        if(pool[clustPop][j,] != edge){
          clusterSize <- clusterSize+1
          clusterMates <- append(clusterMates, as.character(pool[j,]$Number))
        } else {
          break
        }
      }
      clusterMates <- sort(clusterMates)
      cells["clusterSize"][i,] <<- clusterSize
      cells["clusterMates"][i,] <<- paste(unlist(clusterMates), collapse = "_")
    }
  }
  for (i in 1:nrow(cells)){
    if(cells$Number[i] %in% unique(home$Number & cells$clusterID[i] == "None")){
      pool <- subset(home, Number != cells$Number[i])
      pool$lScore <- 
      pCluster <- subset(pool, lScore > clusterCutoff)
      cells$clusterID[cells$Number %in% pCluster$Number,] <- clusterID
      clusterID <- clusterID + 1
    }
  }
}

#-----------------------------------------------------------------
#-----------------------------------------------------------------
#Color palettes:
Color_blind<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cb_black <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#The cellCycle_colors palette are the hexdecimal codes for: 
#                        G1,         S,       G2/M,    mitosis, & premature mitosis respecitvely, as used in Justice et al, 2019
cellCycle_colors <- c("#D4D4D4", "#98C84C", "#23B8CC", "#F16B1A", "#E5001C")



#-----------------------------------------------------------------------------------------------------------
#Program opens required packages
#If these packages aren't install, use the '0_setup.R' function to install them.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(RecordLinkage))
suppressPackageStartupMessages(library(ggbiplot))