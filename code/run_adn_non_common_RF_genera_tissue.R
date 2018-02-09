### Run Random Forest Analysis -- adenoma tissue ALL genera
### Generate model and then test on remaining studies
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "caret", "pROC"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("lu")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer")

# Read in specific data tables to be used

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0)) %>% 
  filter(study %in% c(tissue_sets, both_sets)) %>% 
  rename(sample_id = group) %>% 
  filter(id != 20) #remove the one matched control sample

tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  filter(study %in% c(tissue_sets, both_sets), disease != "cancer") %>% 
  rename(sample_id = group)


####################################################################################################
########################## Functions to run the analysis ###########################################
####################################################################################################











