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





##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


shared_list <- mapply(get_file, c("wang", "brim"), "data/process/", "shared")

taxa_list <- mapply(get_file, c("wang", "brim"), "data/process/", "taxonomy", SIMPLIFY = F)


genera_data <- mapply(get_tax_level_shared, c("wang", "brim"), 
                      "shared_list", "taxa_list", 6)














