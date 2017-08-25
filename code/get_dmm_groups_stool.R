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


# Function to get a subsamples genus file
get_genera_subsample <- function(i, dataList = genera_files){
  
  tempData <- as.matrix(dataList[[i]])
  
  lowest_seq_count <- min(rowSums(tempData))
  
  total_genera <- length(colnames(tempData))
  
  genera_names <- colnames(tempData)
  
  stored_draws_List <- NULL
  
  for(j in rownames(tempData)){
    
    tempdraw <- c()
    tempVector <- unname(tempData[j, ])

    for(k in 1:total_genera){
      
      tempdraw <- c(tempdraw, rep(k, tempVector[k]))
      
    }
    
    stored_draws_List[[j]] <- tempdraw
  }
  
  stored_rd <- lapply(stored_draws_List, function(x) sample(x, lowest_seq_count))
  
  num_counts <- lapply(stored_rd, function(x) as.data.frame(table(x), stringsAsFactors = FALSE))
  
  agg_genera <- lapply(num_counts, 
                 function(x) assign_genera(x, genera_names))
  
  final_table <- t(as.data.frame.list(agg_genera))
  
  return(final_table)
}


# Function to get the assignments needed from the sampling
assign_genera <- function(dataTable, generaVector){
  
  tempVector <- rep(0, length(generaVector))
  names(tempVector) <- generaVector
  
  
  for(i in 1:length(dataTable[, "x"])){
    
    tempVector[as.numeric(dataTable[i, "x"])] <- dataTable[i, "Freq"]
    
  }
  
  
  return(tempVector)  
}




# grabs the counts by row (write as for loop first)
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
  

mapply(get_genera_subsample, c("wang", "brim"))



