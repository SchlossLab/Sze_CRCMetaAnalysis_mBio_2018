### Generate DMM groups
### Use the genera files to create DMM groups
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "DirichletMultinomial"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("wang", "brim", "weir", "ahn", "zeller", "baxter")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")


# Function to read in needed genera file
get_file <- function(i, path_to_file, ending){
  
  temp_shared <- read.csv(paste(path_to_file, i, "/", i, ending, sep = ""), 
                            header = T, stringsAsFactors = F)
  print(paste("Completed Reading in: ", i, ending, "data", sep = ""))
  
  return(temp_shared)
}




##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################



test <- get_file("wang", "data/process/", "_genera_shared.csv")

