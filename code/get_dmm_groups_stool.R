### Generate DMM groups
### Use the genera files to create DMM groups
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "DirichletMultinomial", "vegan"))

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
                            header = T, stringsAsFactors = F, row.names = 1)
  print(paste("Completed Reading in: ", i, ending, " data", sep = ""))
  
  return(temp_shared)
}


# Function to get 
get_genera_subsample <- function(i, dataList = genera_files){
  
  tempData <- dataList[[i]]
  
  lowest_seq_count <- min(rowSums(tempData))
  
  
  return(lowest_seq_count)
}


#lowest_total_seq <- min(rowSums(wang_genera))
# gets the lowest total sequences
# create function that gets the length of columns
# stores column names
# grabs the counts by row
# creates a new vector based on these parameters
# randomly samples this new vector 
    # adds the values up for each one
    # Add option for how many times 
    # get the average 
    # probably needs to insert 0's between columns
# saves this output





### Need to create random sampling function that averages x number of subsamplings


##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################



genera_files <- mapply(get_file, c("wang", "brim"), "data/process/", "_genera_shared.csv")
  





