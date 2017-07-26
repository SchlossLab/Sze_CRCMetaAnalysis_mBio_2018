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


# Function to get matched data and separate them from other data sets

### BURNS
test <- tissue_metadata[["burns"]] %>% select(Run_s, host_subject_id_s) %>% 
  rename(group = Run_s)

alpha_test <- pwr_transformed_data[["burns"]]

combined_data <- test %>% inner_join(alpha_test, by = "group") %>% 
  group_by(host_subject_id_s) %>% filter(n()>1)

non_combined_data <- test %>% inner_join(alpha_test, by = "group") %>% 
  group_by(host_subject_id_s) %>% filter(n() ==1)

### DEJEA
d_test <- tissue_metadata[["dejea"]] %>% select(Run_s, Sample_Name_s) %>% 
  rename(group = Run_s) %>% 
  separate(Sample_Name_s, c("id", "ext_disease", "ext_location"), sep = "\\.")

d_alpha_test <- pwr_transformed_data[["dejea"]]


d_combined_data <- d_test %>% inner_join(d_alpha_test, by = "group") %>% 
  group_by(id) %>% filter(n()>1)

d_non_combined_data <- d_test %>% inner_join(d_alpha_test, by = "group") %>% 
  group_by(id) %>% filter(n() == 1)

### LU









# Read in transformed data
pwr_transformed_data <- mapply(get_pwr_transformed_data, c(tissue_sets, both_sets), 
                               SIMPLIFY = F)

# Read in needed metadata for tissue
tissue_metadata <- mapply(get_orig_metadata, c(tissue_sets, both_sets), SIMPLIFY = F)



# Chen is unmatched - chen.metadata
# Flemer is unmatched - flemer.metadata
# Sanaparedy is unmatched - sana.metadata


# Burns is matched (want host_subject_id_s column) - burnsMetadata.csv
# Lu data A matches with B - luData.csv
# Dejea is matched (host_tissue_sampled_s column) - SraRunTable.txt
# Geng is matched (Sample_Name_s column) - GengData.txt





