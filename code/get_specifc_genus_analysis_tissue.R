### Pull, transform, and normalize 4 crc genera 
### Specifically analyze the genera tied to crc from previous research (tissue)
### One for unmatched and one for matched if possible
### Marc Sze



# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "epiR", "metafor"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("dejea", "geng", "sana", "burns")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")

# CRC genera of interest
crc_genera <- c("Fusobacterium", "Peptostreptococcus", "Porphyromonas", "Parvimonas")



# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_data <- function(i, metadata_table){
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
  study_meta <- metadata_table %>% filter(study == i) %>% rename(sampleID = group)
    
  # conditional that checks for whether length of rows of meta data is smaller
  if(length(rownames(study_meta)) < length(rownames(sub_genera_data))){
    # grab only the samples in the meta data file for down stream analysis
    sub_genera_data <- sub_genera_data %>% slice(match(study_meta$sampleID, sample_ID))
    
  } else{
    # grab only files in the data file for analysis
    study_meta <- study_meta %>% slice(match(sub_genera_data$sample_ID, sampleID))
  }
  # re assigns the rown names while removing the extra column used for sorting
  sample_names <- sub_genera_data$sample_ID
  # creates a list file with both data sets
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta)
  # returns the combined list file
  return(dataList)
  
}



# Function to grab only the genera file and pull specific genus from it
get_specific_genera <- function(i, genera_to_get, table_name, meta_name, 
                                dataList){
  # i represents the study
  # genera_to_get represents the genera of interest to pull specifically
  # table_name represents the table the has all the genera data (subsampled)
  # meta_name represents the metadata name that stores the relevent meta data
  # dataList is a list with both meta data and sub sampled genus data
  
  # grab the specific genera and merge with the meta data file
  tempData <- dataList[[i]][[table_name]] %>% 
    select(sample_ID, one_of(genera_to_get)) %>% 
    rename(sampleID = sample_ID) %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    inner_join((dataList[[i]][[meta_name]] %>% 
                  mutate(sampleID = as.character(sampleID), 
                         disease = ifelse(disease == "normal", 
                                          invisible("control"), invisible(disease)))), 
               by = "sampleID") %>% 
    mutate(disease2 = ifelse(disease != "cancer", invisible("control"), invisible(disease))) %>% 
    as.data.frame()
  
  # Print output when merged completed
  print(paste("Finished pulling and merging files for", i, "data sets"))
  # Return the newly transformed data
  return(tempData)
  
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Read in specific data tables to be used

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0)) %>% 
  select(-sobs, -shannon, -shannoneven, -r_sobs, -r_shannon, -r_shannoneven)

tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  mutate(disease = gsub("adenoma", "polyp", disease)) %>% 
  select(-sobs, -shannon, -shannoneven, -r_sobs, -r_shannon, -r_shannoneven)


# Remove polyp only group
no_p_tissue_matched <- tissue_matched %>% filter(study != "lu")
no_p_tissue_unmatched <- tissue_unmatched %>% filter(study != "lu")

# Generate RR and data tables for every study
ind_matched_data <- sapply(c("burns", "dejea", "geng"), 
                           function(x) get_data(x, tissue_matched), simplify = F)

ind_unmatched_data <- sapply(c("burns", "dejea", "sana", both_sets),  
                             function(x) get_data(x, tissue_matched), simplify = F)


# pull the specific genera of interest and merge with the meta data
matched_specific_genera_list <- sapply(
  c("burns", "dejea", "geng"), 
  function(x) get_specific_genera(x, crc_genera, 
                                  "sub_genera_data", "study_meta", 
                                  ind_matched_data), simplify = F)

unmatched_specific_genera_list <- sapply(
  c("burns", "dejea", "sana", both_sets), 
  function(x) get_specific_genera(x, crc_genera, 
                                  "sub_genera_data", "study_meta", 
                                  ind_unmatched_data), simplify = F)











