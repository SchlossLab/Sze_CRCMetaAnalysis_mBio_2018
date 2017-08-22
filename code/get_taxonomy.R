### Generate genera files
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr"))

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


# Function to read in needed shared or taxonomy file
get_file <- function(i, path_to_file, ending){
  
  temp_shared <- read.delim(paste(path_to_file, i, "/", i, ".", ending, sep = ""), 
                            header = T, stringsAsFactors = F)
  
  return(temp_shared)
}


# Function to remove select columns
get_column <- function(i, dataList_name, keepMethod, keepGroup){
  
  dataList <- get(dataList_name)
  
  tempData <- dataList[[i]] %>% filter(method == keepMethod) %>% 
    select(contains(keepGroup))
  
  return(tempData[, keepGroup])
}


# Function to subset genera files
remake_genera <- function(i, dataList_name, analyzed_name){
  
  dataList <- get(dataList_name)
  analyzedList <- get(analyzed_name)
  
  tempData <- dataList[[i]]
  tempVector <- analyzedList[[i]]
  
  tempData <- tempData %>% slice(match(tempVector, Group))
  
  print(paste("The order is the same:", identical(tempData$Group, tempVector)))
  
  return(tempData)
}


#### Need to do
# Need to limit shared files by sub sampled data that had alpha analyzed 
# save data files to be used later for DMM and other downstream analysis


##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Get shared data
shared_list <- mapply(get_file, c("wang", "brim"), "data/process/", "shared")

# Get taxa data
taxa_list <- mapply(get_file, c("wang", "brim"), "data/process/", "taxonomy", SIMPLIFY = F)

# Get alpha data to get correct groups analyzed
alpha_analyzed_list <- mapply(get_file, c("wang", "brim"), "data/process/", 
                              "groups.ave-std.summary", SIMPLIFY = F)

# Get genera data
genera_data <- mapply(get_tax_level_shared, c("wang", "brim"), 
                      "shared_list", "taxa_list", 6)

# Get analyzed samples by study
analyzed_samples <- mapply(get_column, c("wang", "brim"), "alpha_analyzed_list", "ave", "group")

# subset the genera data by the samples analyzed for alpha
reordered_genera <- mapply(remake_genera, c("wang", "brim"), "genera_data", "analyzed_samples")









