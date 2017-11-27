### Run Random Forest Analysis OTU Level - Tissue
### Generate model and then test on remaining studies
### Marc Sze

## For tissue set CV to 5 since there are in general less samples per study

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "stringr"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
# Ignore weir since none have nzv
stool_sets <- c("wang", "ahn", "zeller", "baxter", "hale")
adn_stool_sets <- c("brim", "zeller", "baxter", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")


##############################################################################################
########################## Group of Functions needed to run the analysis #####################
##############################################################################################


# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_data <- function(i, file_path, ending, genera = F){
  # i is the study of interest
  
  if(genera == F){
    
    temp_taxa <- read_tsv(paste("data/process/", i, "/", i, ".taxonomy", sep = ""))
    # grabs subsampled data and assigns rownames from sample names to table
    tempData <- read_csv(paste(file_path, i, ending, 
                               sep = "")) %>% 
      left_join(temp_taxa, by = c("otu" = "OTU")) %>% 
      mutate(Taxonomy = str_replace_all(Taxonomy, "\\((\\d{2,3})\\)", "")) %>% 
      separate(Taxonomy, c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>% 
      select(-kingdom, -class, -order, -species, -Size) %>% 
      mutate_at(vars(phylum:genus), function(x) str_replace(x, "_unclassified", ""))
    
    
  } else{
    
    # grabs subsampled data and assigns rownames from sample names to table
    tempData <- read_csv(paste(file_path, i, ending, sep = ""))
    
  }
  
  
  # returns the combined list file
  return(tempData)
  
}


# Function to collect the top X OTUs or genera and see how many times it occurs across studies
get_occurances <- function(study_vector, dataList, lowest_var, var_of_int){
  
  tempData <- sapply(study_vector, function(x) slice(dataList[[x]], 1:lowest_var) %>% 
                       mutate(study = x) %>% distinct_(var_of_int, .keep_all = T), simplify = F) %>% 
    bind_rows() %>% 
    select(one_of(var_of_int)) %>% 
    mutate_(var_of_int = str_replace(var_of_int, "Ruminococcus2", "Ruminococcus")) %>% 
    count_(var_of_int) %>% 
    arrange(desc(n)) %>% 
    rename(occurance = n)
  
  return(tempData)
}



##############################################################################################
########################## Code used to run the analysis  ####################################
##############################################################################################

# Read in the summary important Genera tables
crc_imp_genera <- sapply(c(stool_sets, "flemer"), 
                         function(x) get_data(x, "data/process/tables/genus_stool_RF_full_", 
                                              "_imp_vars.csv", genera = T), simplify = F)

adn_imp_genera <- sapply(adn_stool_sets, 
                         function(x) get_data(x, "data/process/tables/adn_genus_stool_RF_full_", 
                                              "_imp_vars.csv", genera = T), simplify = F)


# Read in the summary important OTU tables
crc_imp_otu_data <- sapply(c(stool_sets, "flemer"), 
                       function(x) get_data(x, "data/process/tables/", 
                                            "_imp_otu_table.csv", genera = F), simplify = F)

adn_imp_otu_data <- sapply(adn_stool_sets, 
                           function(x) get_data(x, "data/process/tables/adn_", 
                                                "_imp_otu_table.csv", genera = F), simplify = F)


test <- get_occurances(c(stool_sets, "flemer"), crc_imp_genera, lowest_var = 10, "otu")


test <- get_similarity(adn_stool_sets, adn_imp_genera, lowest_var = 10)
