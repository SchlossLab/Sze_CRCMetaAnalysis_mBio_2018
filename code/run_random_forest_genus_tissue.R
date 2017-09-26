### Run Random Forest Analysis
### Generate model and then test on remaining studies
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "caret", "pROC"))

# Read in specific data tables to be used

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0)) %>% 
  # remove polyp sample
  filter(id != "3776") %>% 
  select(group, id, disease, study) %>% 
  rename(sample_id = group)
  

  
tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  mutate(disease = gsub("adenoma", "polyp", disease)) %>% 
  select(group, id, disease, study) %>% 
  rename(sample_id = group)


# Get matched studies
matched_studies <- unique(tissue_matched$study)
unmatched_studies <- unique(tissue_unmatched$study)





##############################################################################################
########################## Group of Functions needed to run the analysis #####################
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
    select(sample_ID, everything())
  # grabs the meta data and transforms polyp to control (polyp/control vs cancer) 
  study_meta <- metadata %>% filter(study == i)
    
  # Looks for Na in the meta data of interest and removes respective samples
  study_meta <- study_meta %>% filter(!is.na(disease))
  
  # conditional that checks for whether length of rows of meta data is smaller
  if(length(rownames(study_meta)) < length(rownames(sub_genera_data))){
    # grab only the samples in the meta data file for down stream analysis
    sub_genera_data <- sub_genera_data %>% slice(match(study_meta$sample_id, sample_ID))
    
  } else{
    # grab only files in the data file for analysis
    study_meta <- study_meta %>% slice(match(sub_genera_data$sample_ID, sample_id))
  }
  # Prints out the total number of genera for that specific study
  print(paste("Total number of columns in", i, "is", 
              length(colnames(sub_genera_data))))
  # re assigns the rown names while removing the extra column used for sorting
  sample_names <- sub_genera_data$sample_ID
  # creates a list file with both data sets
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













##############################################################################################
########################## Code used to run the analysis (unmatched) #########################
##############################################################################################


# reads in all the stool data into one list
unmatched_stool_study_data <- sapply(unmatched_studies, 
                           function(x) get_data(x, tissue_unmatched), simplify = F)


# Align the genera so there is the same number for each data set.
unmatched_matched_genera_list <- align_genera(unmatched_studies, "column_length", 
                                    "sub_genera_data", unmatched_stool_study_data)









