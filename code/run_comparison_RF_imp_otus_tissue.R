### GEnerate Summary data for important genus and OTU - Tissue
### For both adenoma and CRC
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "stringr"))

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

# Get studies that contain matched samples
# Remove Lu since it only has polyps
matched_studies <- unique(
  tissue_matched$study[!(tissue_matched$study %in% c("lu"))])

# Get studies that contain unmatched samples
# Need to remove dejea and lu
# Lu only has polyps
# Dejea only has cancer 
unmatched_studies <- unique(
  tissue_unmatched$study[!(tissue_unmatched$study %in% c("dejea", "lu"))]) 

adn_tissue <- c("flemer", "lu")

rm(tissue_matched, tissue_unmatched)

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
      mutate_at(vars(phylum:genus), function(x) str_replace(x, "(\\d$)", "")) %>% 
      mutate(genus = ifelse(grepl("", genus) == T, invisible(family), invisible(genus)))
    
    
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
crc_unmatched_imp_genera <- sapply(unmatched_studies, 
                         function(x) get_data(x, "data/process/tables/ALL_genus_unmatched_tissue_RF_", 
                                              "_imp_vars.csv", genera = T), simplify = F)

crc_matched_imp_genera <- sapply(matched_studies, 
                         function(x) get_data(x, "data/process/tables/ALL_genus_matched_tissue_RF_", 
                                              "_imp_vars.csv", genera = T), simplify = F)


adn_imp_genera <- sapply(adn_tissue, 
                                 function(x) get_data(x, "data/process/tables/adn_ALL_genus_unmatched_tissue_RF_full_", 
                                                      "_imp_vars.csv", genera = T), simplify = F)

# Read in the summary important OTU tables
crc_unmatched_imp_otu_data <- sapply(unmatched_studies, 
                           function(x) get_data(x, "data/process/tables/", 
                                                "_unmatched_tissue_imp_otu_table.csv", genera = F), simplify = F)

crc_matched_imp_otu_data <- sapply(matched_studies, 
                           function(x) get_data(x, "data/process/tables/", 
                                                "_matched_tissue_imp_otu_table.csv", genera = F), simplify = F)

adn_imp_otu_data <- sapply(adn_tissue, 
                                   function(x) get_data(x, "data/process/tables/adn_", 
                                                        "_tissue_imp_otu_table.csv", genera = F), simplify = F)


# Generate the occurances tables
crc_unmatched_genera_occurances <- get_occurances(unmatched_studies, crc_unmatched_imp_genera, lowest_var = 10, "otu")
crc_matched_genera_occurances <- get_occurances(matched_studies, crc_matched_imp_genera, lowest_var = 10, "otu")
adn_genera_occurances <- get_occurances(adn_tissue, adn_imp_genera, lowest_var = 10, "otu")

crc_unmatched_otu_occurances <- get_occurances(unmatched_studies, crc_unmatched_imp_otu_data, lowest_var = 10, "genus")
crc_matched_otu_occurances <- get_occurances(matched_studies, crc_matched_imp_otu_data, lowest_var = 10, "genus")
adn_otu_occurances <- get_occurances(adn_tissue, adn_imp_otu_data, lowest_var = 10, "genus")

# Read out data tables
write_csv(crc_unmatched_genera_occurances, "data/process/tables/crc_RF_genera_unmatched_tissue_top10.csv")
write_csv(crc_matched_genera_occurances, "data/process/tables/adn_RF_genera_matched_tissue_top10.csv")
write_csv(adn_genera_occurances, "data/process/tables/adn_RF_genera_tissue_top10.csv")
write_csv(crc_unmatched_otu_occurances, "data/process/tables/crc_RF_otu_unmatched_tissue_top10.csv")
write_csv(crc_matched_otu_occurances, "data/process/tables/adn_RF_otu_matched_tissue_top10.csv")
write_csv(adn_otu_occurances, "data/process/tables/adn_RF_otu_tissue_top10.csv")
