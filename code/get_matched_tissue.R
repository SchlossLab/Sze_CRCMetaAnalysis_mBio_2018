### Get Matched Tissue samples
### For tissue comparisons only want to keep matched samples
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "rcompanion"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")

# Create vector with data that has matched samples
matched_sets <- c("burns", "dejea", "lu", "geng")



# Function to read in each power transformed data sets
get_pwr_transformed_data <- function(i){
  
  data_table <- read.csv(
    paste("data/process/tables/", i, "_tissue_transformed_data_alpha_raw_values.csv", sep = ""), 
    header = T, stringsAsFactors = F)
  
  return(data_table)
}


# Function to read in complete tissue metadata
get_orig_metadata <- function(i){
  
  if(i == "chen" | i == "flemer" | i == "sana"){
    
    data_table <- read.delim(
      paste("data/process/", i, "/", i, ".metadata", sep = ""), 
      header = T, stringsAsFactors = F)
  } else if(i == "dejea" | i == "geng" | i == "lu"){
    
    if(i == "dejea" | i == "geng"){
      
      download_name <- list.files(path = paste("data/process/", i, sep = ""), 
                                  pattern = "*.txt")
      
      data_table <- read.delim(paste("data/process/", i, "/", download_name, sep = ""), 
                               header = T, stringsAsFactors = F)
    } else {
      
      download_name <- list.files(path = paste("data/process/", i, sep = ""), 
                                  pattern = "*.csv")
      
      data_table <- read.csv(paste("data/process/", i, "/", download_name, sep = ""), 
                               header = T, stringsAsFactors = F)
    } 
    
  } else {
    
    download_names <- list.files(path = paste("data/process/", i, sep = ""), 
                                pattern = "*.csv")
    wanted_download <- unique(gsub("[0-9]", "", download_names))
    
    data_table <- read.csv(paste("data/process/", i, "/", wanted_download, sep = ""), 
                           header = T, stringsAsFactors = F)
  }
  
  return(data_table)
}


# Function to apply a z-score normalization and and a new column with these values
# for sobs, shannon, and shannoneven
zscore_transform <- function(dataList, 
                             variables = c("sobs", "shannon", "shannoneven")){
  # dataList represents the list data that was read in by get_transformed_data
  # variables reprents data of interest for downstream comparisons
  
  # command that z-score transforms interested data 
  z_trans_data <- dataList %>% 
    mutate_at(vars(variables), funs(as.numeric(scale(.))))
  
  # return to working environment the new z-score transformed data
  return(z_trans_data)
}



# Function to get matched data and separate them from other data sets
get_combinations <- function(i, matching = data_to_match_list, 
                             alpha_set = alpha_to_match_list){
  data_list <- c()
  
  data_list[["matched"]] <- matching[[i]] %>% 
    inner_join(alpha_set[[i]], by = "group") %>% 
    group_by(id) %>% filter(n() > 1)
  
  data_list[["unmatched"]] <- matching[[i]] %>% 
    inner_join(alpha_set[[i]], by = "group") %>% 
    group_by(id) %>% filter(n() == 1)
  
  return(data_list)
}


# combine the non-matched data with other data sets and create a stand alone matched data set




# Read in transformed data
pwr_transformed_data <- mapply(get_pwr_transformed_data, c(tissue_sets, both_sets), 
                               SIMPLIFY = F)

# Z-Score normalize the data
zscore_pwr_transform_data <- lapply(pwr_transformed_data, 
                                    function(x) zscore_transform(x))



# Read in needed metadata for tissue
tissue_metadata <- mapply(get_orig_metadata, c(tissue_sets, both_sets), SIMPLIFY = F)


# Create needed lists of data
data_to_match_list <- list(
  burns = tissue_metadata[["burns"]] %>% select(Run_s, host_subject_id_s) %>% 
    rename(id = host_subject_id_s, group = Run_s), 
  dejea = tissue_metadata[["dejea"]] %>% select(Run_s, Sample_Name_s) %>% 
    rename(group = Run_s) %>% 
    separate(Sample_Name_s, c("id", "ext_disease", "ext_location"), sep = "\\."), 
  lu = tissue_metadata[["lu"]] %>% 
    mutate(Sample_Name_s = gsub("[A-Z]", "", Sample_Name_s)) %>% 
    rename(id = Sample_Name_s, group = Run_s), 
  geng = tissue_metadata[["geng"]] %>% 
    select(Run_s, Sample_Name_s) %>% 
    mutate(Sample_Name_s = gsub("[A-Z]", "", Sample_Name_s)) %>% 
    rename(id = Sample_Name_s, group = Run_s))

alpha_to_match_list <- zscore_pwr_transform_data[matched_sets]

# generate the combined data
combination_data <- mapply(get_combinations, matched_sets, USE.NAMES = T, SIMPLIFY = F)





# Chen is unmatched - chen.metadata
# Flemer is unmatched - flemer.metadata
# Sanaparedy is unmatched - sana.metadata


# Burns is matched (want host_subject_id_s column) - burnsMetadata.csv
# Lu data A matches with B - luData.csv
# Dejea is matched (host_tissue_sampled_s column) - SraRunTable.txt
# Geng is matched (Sample_Name_s column) - GengData.txt





