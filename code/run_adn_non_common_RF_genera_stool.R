### Run Random Forest Analysis -- Adenoma
### Generate model and then test on remaining studies
### This code manually inserts data sets after the nzv step
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "caret", "pROC"))


# Stool Only polyp sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("brim", "zeller", "baxter", "hale")


####################################################################################################
########################## Functions to run the analysis ###########################################
####################################################################################################


# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_data <- function(i){
  # i is the study of interest
  
  # gets original sample names
  sample_names <- (get_file(i, "data/process/", "_genera_shared.csv") %>% 
                     mutate(sample_names = rownames(.)) %>% select(sample_names))[, "sample_names"]
  # grabs subsampled data and assigns rownames from sample names to table
  sub_genera_data <- get_file(i, "data/process/", "_subsample_genera.csv", rows_present = F, 
                              vec_of_rownames = sample_names) %>% 
    as.data.frame() %>% mutate(sample_ID = rownames(.)) %>% 
    select(sample_ID, everything()) %>% mutate(sample_ID = as.character(sample_ID))
  # grabs the meta data and transforms polyp to control (polyp/control vs cancer) 
  study_meta <- get_file(i, "data/process/", ".metadata", rows_present = F,  
                         "stool", metadata = T) %>% 
    filter(disease != "cancer", !is.na(disease)) %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    select(sampleID, disease)
  
  sub_genera_data <- study_meta %>% 
    inner_join(sub_genera_data, by = c("sampleID" = "sample_ID")) %>% 
    select(-sampleID)
  
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta, 
                   column_length = length(colnames(sub_genera_data)))
  # returns the combined list file
  return(dataList)
  
}


####################################################################################################
########################## Functions to run the analysis ###########################################
####################################################################################################

# reads in all the stool data into one list
stool_study_data <- mapply(get_data, stool_sets, SIMPLIFY = F)





