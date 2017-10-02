### Run Random Forest Analysis -- adenoma tissue
### Generate model and then test on remaining studies
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "caret", "pROC"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("lu")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer")

# Read in specific data tables to be used

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0)) %>% 
  filter(study %in% c(tissue_sets, both_sets)) %>% 
  rename(sample_id = group) %>% 
  filter(id != 20) #remove the one matched control sample

tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  filter(study %in% c(tissue_sets, both_sets), disease != "cancer") %>% 
  rename(sample_id = group)


##############################################################################################
############### List of functions to run the analysis ########################################
##############################################################################################

# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_data <- function(i, metadata){
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
  study_meta <- metadata %>% filter(study == i)
  
  # Looks for Na in the meta data of interest and removes respective samples
  study_meta <- study_meta %>% filter(!is.na(disease)) %>% 
    select(sample_id, disease)
  
  sub_genera_data <- study_meta %>% 
    inner_join(sub_genera_data, by = c("sample_id" = "sample_ID")) %>% 
    select(-sample_id)
  
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta, 
                   column_length = length(colnames(sub_genera_data)))
  # returns the combined list file
  return(dataList)
  
}


# Function that aligns the genera from all data sets
align_genera <- function(studies, length_column_name, 
                         genera_data_name, dataList){
  # studies is a variable for the name of study to be used
  # length_column_name is the name of the vector in dataList that stores the total genera number
  # genera_data_name is the name of the data table in dataList that contains the genera information
  # dataList is a list that has for every study the genus info, metadata, and total genera present
  
  # Pulls out the total number of genera identified in each study and orders them highest to lowest
  genera_num_list <- sort.int(sapply(studies, 
                                     function(x) dataList[[x]][[length_column_name]]))
  # Counting variable (helps to direct flow)
  x = 1
  # Checks to see if the lowest value is the same as the highest, if not keep iterating through
  while(genera_num_list[1]!= genera_num_list[length(genera_num_list)]){
    # stores the  lowest genera study name
    lowest_genera_study <- names(genera_num_list[1])
    # Check to see if this is the first time through the iteration
    if(x == 1){
      # Remove any groups that have unclassified as part of their ID
      genera_names <- dataList[[lowest_genera_study]][[genera_data_name]] %>% 
        select(-contains("_unclassified")) %>% colnames(.)
    } else{
      # If not the first through get the names of genrea in the current lowest study 
      genera_names <- colnames(temp_aligned_genera[[lowest_genera_study]])
    }
    
    # iterate through each study matching only those in the lowest genera data set
    temp_aligned_genera <- suppressWarnings(sapply(studies, 
                                                   function(x) 
                                                     select(dataList[[x]][[genera_data_name]], 
                                                            one_of(genera_names)), simplify = F))
    # get the updated total genera and then sort lowest to highest
    genera_num_list <- sort.int(sapply(studies, 
                                       function(x) 
                                         length(colnames(temp_aligned_genera[[x]]))))
    # Print out an update as to how the matching is going
    print(paste("Min and Max total genera is:", 
                min(genera_num_list), ",", max(genera_num_list)))
    # move the tracker / flow director up one
    x = x + 1
    
  }
  # When the while loop exits return the aligned genera data list
  return(temp_aligned_genera)
}


# Function that grabs the meta data and replaces sampleID with disease call
assign_disease <- function(studies, matched_genera){
  # studies is the variable with the names of the studies 
  # matched_genera is the data list with only genera in every study
  
  # Get the respective metadata file of interest
  tempData <- matched_genera[[studies]]
  # Gets transforms sample_ID column into a disease column with control v cancer calls
  corr_tempData <- tempData %>% 
    mutate(disease = factor(ifelse(disease == "normal", 
                                   invisible("control"), invisible(disease)), 
                            levels = c("control", "polyp")))
  # Returns the modified data frame that can be used for RF analysis
  return(corr_tempData)
  
}




##############################################################################################
########################## Code used to run the analysis (unmatched) #########################
##############################################################################################

# reads in all the stool data into one list
unmatched_stool_study_data <- sapply(c(both_sets, tissue_sets), 
                                     function(x) get_data(x, tissue_unmatched), simplify = F)

#Align the genera so there is the same number for each data set.
unmatched_matched_genera_list <- align_genera(c(both_sets, tissue_sets), "column_length", 
                                    "sub_genera_data", unmatched_stool_study_data)

# Generate data sets to be used in random forest
unmatched_rf_datasets <- sapply(c(both_sets, tissue_sets), 
                      function(x) assign_disease(x, unmatched_matched_genera_list), simplify = F)


##############################################################################################
########################## Code used to run the analysis (unmatched) #########################
##############################################################################################

# reads in all the stool data into one list
matched_stool_study_data <- sapply(c(tissue_sets), 
                                   function(x) get_data(x, tissue_matched), simplify = F)






