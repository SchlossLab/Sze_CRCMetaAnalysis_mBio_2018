### Get Tissue relative risk comparisons
### One for unmatched and one for matched if possible
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "epiR", "metafor"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")




# Function to read in the respective transformed data tables for each study
# These should have been previously power transformed to correct skew
get_data <- function(i, sampleType){
  # i represents the data set
  # sampleType represents whether it is stool or tissue
  
  # Command that actually does the reading in
  data_list <- read.csv(paste("data/process/tables/", i, "_", sampleType, "_",  
                              "alpha_raw_values.csv", sep = ""), 
                        header = T, stringsAsFactors = F)
  
  # return to working environment the data list
  return(data_list)
}








# Read in the respective data
tissue_data <- mapply(get_data, c(tissue_sets, both_sets), "tissue", SIMPLIFY = F)








