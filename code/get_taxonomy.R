### Generate genera files
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "foreach", "doMC"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# ignore Hale and Zeller for testing since the files are too large
stool_sets <- c("wang", "brim", "weir", "ahn","zeller", "baxter")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore flemer since it is too large right now
both_sets <- c("flemer", "chen")


# Function to read in needed shared or taxonomy file
get_file <- function(i, path_to_file, ending){
  
  temp_shared <- read.delim(paste(path_to_file, i, "/", i, ".", ending, sep = ""), 
                            header = T, stringsAsFactors = F)
  print(paste("Completed Reading in:", i, ending, "data"))
  
  return(temp_shared)
}


# Function to remove select columns
get_column <- function(dataList_name, keepMethod, keepGroup){
  
  dataList <- dataList_name
  
  tempData <- dataList %>% filter(method == keepMethod) %>% 
    select(contains(keepGroup))
  
  return(tempData[, keepGroup])
}


# Function to subset genera files
remake_genera <- function(dataList, analyzedList){
  
  tempData <- dataList
  tempVector <- analyzedList
  
  tempData <- tempData %>% slice(match(tempVector, Group))
  
  print(paste("The order is the same:", identical(tempData$Group, tempVector)))
  
  return(tempData)
}


# Function to write out the new genera data

save_data <- function(i, genera_file, path_to_file, ending){
  
  write.csv(genera_file, 
            paste(path_to_file, i, "/", i, ending, sep = ""), 
            row.names = F)
  
  print(paste("Completed Saving", i, "taxa data"))
}


# Master Control Function
run_get_genera <- function(i){
  
  shared <- get_file(i, "data/process/", "shared")
  taxa <- get_file(i, "data/process/", "taxonomy")
  alpha_data <- get_file(i, "data/process/", "groups.ave-std.summary")
  
  genera_data <- get_tax_level_shared(i, shared, taxa, 6)
  analyed_data <- get_column(alpha_data, "ave", "group")
  reordered_genera <- remake_genera(genera_data, analyed_data)
  
  save_data(i, reordered_genera, "data/process/", "_genera_shared.csv")
}




##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Call number of processors to use
#registerDoMC(cores = 8)

for(i in c(stool_sets, tissue_sets, both_sets)){

  run_get_genera(i)
}

#foreach(i = c(stool_sets, tissue_sets, both_sets)) %dopar% run_get_genera(i)





