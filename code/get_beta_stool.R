### Get Stool Beta Diversity (Bray-Curtis)
### Look for overall group differences
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "lme4"))

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("wang", "brim", "weir", "ahn", "zeller", "baxter")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")


# Create function to read in distance data
get_distance <- function(i, path_to_file, file_ending){
  
  tempDist <- read.dist(paste(path_to_file, i, "/", i, ".", file_ending, sep = ""))
  
  return(tempDist)
}


# This function below is used as a proxy of the full metadata table since only 
# need a small component and it is contained in these files
get_metadata_data <- function(i, sampleType){
  # i represents the data set
  # sampleType represents whether it is stool or tissue
  
  # Command that actually does the reading in
  data_list <- read.csv(paste("data/process/tables/", i, "_", sampleType, "_",  
                              "transformed_data_alpha_raw_values.csv", sep = ""), 
                        header = T, stringsAsFactors = F) %>% 
    select(-sobs, -shannon, -shannoneven)
  
  # return to working environment the data list
  return(data_list)
}






# Read in distance data
distance_data <- mapply(get_distance, c(stool_sets, both_sets), 
                        "data/process/", "braycurtis.0.03.lt.ave.dist")

# Read in data with needed metadata
metadata <- mapply(get_metadata_data, c(stool_sets, both_sets), "stool", SIMPLIFY = F)


#### Need to do list
###### load data in
###### Use distance matrices to generate PERMANOVA values
###### Think of potential way to pool this information together
    ##### E.g. distance of centroids from each other...





