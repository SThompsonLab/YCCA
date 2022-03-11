defaultW <- getOption("warn")
options(warn = -1)
partial_amount <- 1000

scripts <- data.frame(
  "Number" = c(0, 1, 2, 3, 4, 5, 6),
  "Script_Name" = c("Setup",
                    "dirGen",
                    "purgo",
                    "reName",
                    "bud_all",
                    "reset names",
                    "analyze")
)
print("Welcome to MicroCyte. Be sure to have edited your schema file before beginning. Here are your options:")
print(scripts)

while (TRUE) {
  option <- readline(prompt = paste0("What would you like to run (0-", max(scripts$Number), ", or 'end'): "))
  if (option == "end"){
    print("Thank you for using MicroCyte.")
    break
  } else if (option == "0"){
    
# Runs set up script
    print("Running setup now...")
    source("bin/scripts/scripts/0_setup.R")
    print("Set up complete")
  } else if (option == "1"){
    
# Runs Dirgen only
    print("Generating directories...")
    source("bin/scripts/scripts/1_dirGen.R")
    print("Directories complete")
  } else if (option == "2"){
    
# Runs Purgo and asks if you want to rename
    print("Beginning the purge...")
    source("bin/scripts/scripts/2_purgo.R")
    print("Purge complete")
    contun <- readline(prompt = "Would you like to now run reName (Y/n)? ")
    if (contun != "n"){
      print("Renaming images per the schema file...")
      source("bin/scripts/scripts/3_reName.R")
      reName()
      print("Images renamed. Please run YeastTakes to continue.")
    }
  } else if (option == "3"){
    
# Runs ReName only
    print("Renaming images per the schema file...")
    source("bin/scripts/scripts/3_reName.R")
    reName()
    print("Images renamed. Please run YeastTakes to continue.")
  } else if (option == "4"){
    
# Runs Imagen and starts a finish chain
    print("Running Bud et al now...")
    source("bin/scripts/scripts/5_budd.R")

    print("Combining image datasets now...")
    bud()
    print("Combining sample datasets now...")
    budCat()
    print("Datasets combined. Sourcing sirmixaplot for graphing.")
    source("bin/scripts/analysis/SirMixaPlot.R")
    
    print("Please gate populations and then run bud_explore()")
    sirmixaplot("data/cells.csv")

  } else if (option == "5"){
    source("bin/scripts/scripts/3_reName.R")
    reset_names()
  } else if (option == "6"){
    scriptList <- list.files("bin/scripts/analysis/", pattern = ".R")
    print("Sourcing the analysis scripts...")
    cat("\n")
    for(scriptin in scriptList){
      cat(paste0("Sourcing script: ", scriptin))
      cat("\n")
      source(paste0("bin/scripts/analysis/", scriptin))
    }
    cat("\n")
    print("Good luck with your analysis!")
    break
  } else {
    print("Sorry, I didn't catch that. Please try again, or type 'end' to finish.")
    print(scripts)
  }
}

options(warn = defaultW)