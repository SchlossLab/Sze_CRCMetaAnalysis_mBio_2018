### Run Random Forest Analysis OTU Level - Tissue
### Generate model and then test on remaining studies
### Marc Sze

## For tissue set CV to 5 since there are in general less samples per study

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "caret", "pROC"))

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




