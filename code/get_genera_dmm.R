### Get the dmm groups 
### Merge the group data with relevant metadata
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "vegan", "DirichletMultinomial"))

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



# Function to read in needed genera file
get_file <- function(i, path_to_file, ending, rows_present=T, name_of_rows=1, 
                     vec_of_rownames = NULL, make_matrix = F){
  
  if(rows_present == T){
    
    temp_shared <- read.csv(paste(path_to_file, i, "/", i, ending, sep = ""), 
                            header = T, stringsAsFactors = F, row.names = name_of_rows)
  } else{
    
    temp_shared <- read.csv(paste(path_to_file, i, "/", i, ending, sep = ""), 
                            header = T, stringsAsFactors = F)
    rownames(temp_shared) <- vec_of_rownames
  }
  
  
  
  print(paste("Completed Reading in: ", i, ending, " data", sep = ""))
  
  return(temp_shared)
}



# Function to grab groupings from DMM
grab_dmm_groups <- function(dataTable, metaData, kvalue = 2, seed_value = 1234567){
  
  dmm_test <- dmn(dataTable, k = kvalue, verbose = T, seed = seed_value)
  
  dmm_groups <- dmm_test@group
  
  
  
}










##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

for(i in c("wang")){
  
  sample_names <- (get_file(i, "data/process/", "_genera_shared.csv") %>% 
    mutate(sample_names = rownames(.)) %>% select(sample_names))[, "sample_names"]
  
  sub_genera_data <- apply(as.matrix(
    get_file(i, "data/process/", "_subsample_genera.csv", rows_present = F, 
             vec_of_rownames = sample_names)), 2, function(x) round(x))
  
  study_meta <- read.delim(paste("data/process/", i, "/", i, ".metadata", sep = ""), 
                           header = T, stringsAsFactors = F)
  
  
}
