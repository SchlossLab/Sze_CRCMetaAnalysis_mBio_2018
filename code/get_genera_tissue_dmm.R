### Get the dmm groups - tissue
### Merge the group data with relevant metadata
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "vegan", "foreach", "doParallel"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
# Manually add geng where needed
tissue_sets <- c("dejea", "sana", "burns")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "flemer")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")



# Function to create a list of sampleIDs and disease
make_list <- function(i, dataTable_name){
  
  dataTable <- get(dataTable_name)
  
  new_table <- dataTable %>% filter(study == i) %>% 
    select(group, id, disease, sample_type)
  
  if(length(rownames(new_table)) != 0){
    
    return(new_table)
  }
  
}


# Function to read in needed genera file
get_file <- function(i, path_to_file, ending, rows_present=T, name_of_rows=1, 
                     vec_of_rownames = NULL, sample_source, metadata = F){
  # i represents the study name
  # path_to_file is the is the directory where the main studies are stored
  # ending is the extension after i of the name of the file
  # rows_present signifies whether rownames are part of the inputed file
  # name_of_rows represents what column the row names are in
  # sample_source signifies whether it is stool or tissue
  # metadata signifies whether the inputed file is metadata or not
  
  # First conditional to make sure only data goes through this part
  if(metadata == F){
    # second conditional to check if rownames are present or not
    if(rows_present == T){
      # reads in data (assumes csv file) with row names from data table
      temp_shared <- read.csv(paste(path_to_file, i, "/", i, ending, sep = ""), 
                              header = T, stringsAsFactors = F, row.names = name_of_rows)
    } else{
      # reads in data (assumes csv file) with row names from an inputed vector 
      temp_shared <- read.csv(paste(path_to_file, i, "/", i, ending, sep = ""), 
                              header = T, stringsAsFactors = F)
      rownames(temp_shared) <- vec_of_rownames
    }
    
  } else{
    # reads in data table (assumes tab delimited) for meta data
    temp_shared <- read.delim(paste(path_to_file, i, "/", i, ending, sep = ""), 
                              header = T, stringsAsFactors = F)
    
    # conditional that checks for column named sample_type (created during mothur processing)
    if(!("sample_type" %in% colnames(temp_shared))){
      # Create a new column called sample_type if it is not already present
      temp_shared <- temp_shared %>% mutate(sampleID = sample) %>% 
        select(-sample) %>% select(sampleID, everything())
      
    } else{
      # reads in data (assumes tab delimited) if sample_type present
      temp_shared <- temp_shared %>% 
        filter(sample_type == "stool") %>% mutate(sampleID = sample) %>% 
        select(-sample) %>% select(sampleID, everything())
    }
    
  }
  # prints message to stdout updating on completion
  print(paste("Completed Reading in: ", i, ending, " data", sep = ""))
  # returns the necessary temp_shared file
  return(temp_shared)
}


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
    select(sample_ID, everything())

  # returns the combined list file
  return(sub_genera_data)
  
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


# Read in specific data tables with samples that are matched and unmatched

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0))

tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  mutate(disease = gsub("adenoma", "polyp", disease))


matched_meta <- mapply(make_list, c("dejea", "burns", "geng"), "tissue_matched", SIMPLIFY = F)
unmatched_meta <- mapply(make_list, c(tissue_sets, both_sets), "tissue_unmatched", SIMPLIFY = F)

matched_data <- mapply(get_data, c("dejea", "burns", "geng"))
unmatched_data <- mapply(get_data, c(tissue_sets, both_sets))

rm(tissue_matched, tissue_unmatched)














