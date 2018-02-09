### GEnerate Summary data for important genus and OTU - Stool
### For both adenoma and CRC
### Marc Sze

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
      mutate_at(vars(phylum:genus), function(x) str_replace(x, "_unclassified", "")) %>% 
      mutate_at(vars(phylum:genus), function(x) str_replace(x, "(\\d$)", ""))
    
    
  } else{
    
    # grabs subsampled data and assigns rownames from sample names to table
    tempData <- read_csv(paste(file_path, i, ending, sep = "")) %>% 
      mutate(otu = str_replace(otu, "(\\d$)", ""))
    
  }
  
  
  # returns the combined list file
  return(tempData)
  
}


# Function to collect the top X OTUs or genera and see how many times it occurs across studies
get_occurances <- function(study_vector, dataList, lowest_var, var_of_int){
  # study_vector is the vector of studies to analyze
  # dataList is a variable with the read in data
  # lowest_var is the maximum number of OTUs/genera to include
  # var_of_int is the column of interest
  
  tempData <- sapply(study_vector, function(x) slice(dataList[[x]], 1:lowest_var) %>% 
                       mutate(study = x) %>% distinct_(var_of_int, .keep_all = T), simplify = F) %>% 
    bind_rows() %>% 
    select(one_of(var_of_int)) %>% 
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
                         function(x) get_data(x, "data/process/tables/ALL_genus_stool_RF_full_", 
                                              "_imp_vars.csv", genera = T), simplify = F)

adn_imp_genera <- sapply(adn_stool_sets, 
                         function(x) get_data(x, "data/process/tables/adn_ALL_genus_stool_RF_full_", 
                                              "_imp_vars.csv", genera = T), simplify = F)


# Read in the summary important OTU tables
crc_imp_otu_data <- sapply(c(stool_sets, "flemer"), 
                       function(x) get_data(x, "data/process/tables/", 
                                            "_imp_otu_table.csv", genera = F), simplify = F)

adn_imp_otu_data <- sapply(adn_stool_sets, 
                           function(x) get_data(x, "data/process/tables/adn_", 
                                                "_imp_otu_table.csv", genera = F), simplify = F)

# Generate the occurances tables
crc_genera_occurances <- get_occurances(c(stool_sets, "flemer"), crc_imp_genera, lowest_var = 10, "otu")
adn_genera_occurances <- get_occurances(adn_stool_sets, adn_imp_genera, lowest_var = 10, "otu")

crc_otu_occurances <- get_occurances(c(stool_sets, "flemer"), crc_imp_otu_data, lowest_var = 10, "genus")
adn_otu_occurances <- get_occurances(adn_stool_sets, adn_imp_otu_data, lowest_var = 10, "genus")

# Read out data tables
write_csv(crc_genera_occurances, "data/process/tables/crc_RF_genera_stool_top10.csv")
write_csv(adn_genera_occurances, "data/process/tables/adn_RF_genera_stool_top10.csv")
write_csv(crc_otu_occurances, "data/process/tables/crc_RF_otu_stool_top10.csv")
write_csv(adn_otu_occurances, "data/process/tables/adn_RF_otu_stool_top10.csv")






