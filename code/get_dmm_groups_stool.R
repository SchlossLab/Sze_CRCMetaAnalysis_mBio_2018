### Generate DMM groups
### Use the genera files to create DMM groups
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "vegan"))

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
get_genera_subsample <- function(dataList){
  
  tempData <- as.matrix(dataList)
  
  lowest_seq_count <- min(rowSums(tempData))
  
  total_genera <- length(colnames(tempData))
  
  genera_names <- colnames(tempData)
  
  stored_draws_List <- NULL
  
  # Iteratres through each sample in data set
  for(j in rownames(tempData)){
    
    tempdraw <- c()
    tempVector <- unname(tempData[j, ])
    # creates a temp vector with each genus (as a number)
    # repeated based on the number of counts in data frame
    for(k in 1:total_genera){
      
      tempdraw <- c(tempdraw, rep(k, tempVector[k]))
      
    }
    # saves the created data in a new list
    stored_draws_List[[j]] <- tempdraw
  }
  
  # Applies a randome sampling accross the store vector of repeats
  stored_rd <- lapply(stored_draws_List, function(x) sample(x, lowest_seq_count))
  # Generates the counts from the sampling
  num_counts <- lapply(stored_rd, function(x) as.data.frame(table(x), stringsAsFactors = FALSE))
  # Runs the assign genera function
  agg_genera <- lapply(num_counts, 
                 function(x) assign_genera(x, genera_names))
  # Converts the output to a matrix of the same orientation as the data files
  # in the inputted data list
  final_table <- t(as.data.frame.list(agg_genera))
  # Returns the data
  return(final_table)
}


# Function to get the assignments needed from the sampling
assign_genera <- function(dataTable, generaVector){
  # dataTable is part of a list where each dataTable is an individual sample
  # generaVector is a vector of genera names
  
  # Creates a temporary variable of all taxa with 0 and names the vector
  tempVector <- rep(0, length(generaVector))
  names(tempVector) <- generaVector
  
  # Iteratres through each genera sampled and changes the 0 to the
  # correct number of counts
  updatedVector <- lapply(c(1:length(dataTable[, "x"])), 
                       function(x) grab_value(x, tempVector, dataTable))
  

  # returns the final vector
  return(updatedVector)  
}


# Function to pull specific value
grab_value <- function(i, vec_of_int, refTable){
  
  vec_of_int[as.numeric(refTable[i, "x"])] <- refTable[i, "Freq"]
  
  return(vec_of_int)
  
}


# Function to run the sampling x number of times and generate an average from this
get_average_counts <- function(i, repeats, dataList){
  
  total_samples <- length(rownames(dataList))
  genera_names <- colnames(dataList)
  
  full_100_runs <- lapply(1:repeats, function(x) get_genera_subsample(dataList))
  
  
  temp_avg_list <- lapply(c(1:total_samples), 
                          function(x) grab_row(full_100_runs, x, i, genera_names))
  

  final_avg_data <- t(as.data.frame.list(temp_avg_list))
  rownames(final_avg_data) <- rownames(dataList)
  
  print(paste("Completed study ", i, ": taxa subsampling.", sep = ""))
  
  return(final_avg_data)
  
}


# Function to grab rows for averaging 
grab_row <- function(list_of_int, j, study, genera_file){
  
  test <- lapply(list_of_int, function(x) x[j, ])
  
  test <- t(as.data.frame.list(test))
  
  colnames(test) <- genera_file
  rownames(test) <- c(1:length(rownames(test)))
  
  average_vector <- colMeans(test)
  
  return(average_vector)
}


# Function to write out the data
make_file <- function(datafile, path_to_file, ending){
  
  write.csv(datafile, paste(path_to_file, i, "/", i, ending, sep = ""), 
                          row.names = F)
  
  print(paste("Completed writing: ", i, ending, " to file", sep = ""))
}


##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

for(i in c("wang", "ahn")){
  
  genera_data <- get_file(i, "data/process/", "_genera_shared.csv")
  
  avg_subsample_table <- get_average_counts(i, 100, genera_data)
  
  make_file(avg_subsample_table, "data/process/", "_subsample_genera.csv")
  
  rm(genera_data, avg_subsample_table)
}












