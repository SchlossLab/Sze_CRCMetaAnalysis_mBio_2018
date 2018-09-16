### Generate total n used for each study
### Use the RR based count to generate this data
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

# Read in adenoma stool files
adn_stool <- read_csv("data/process/tables/alpha_adn_group_counts_summary.csv") %>% 
  filter(measure == "shannon") %>% 
  mutate(control = high_N + low_N, adenoma = high_Y + low_Y) %>% 
  select(study, control, adenoma)

# Read in the carcinoma stool files
crc_stool <- read_csv("data/process/tables/alpha_group_counts_summary.csv") %>% 
  filter(measure == "shannon") %>% 
  mutate(control = high_N + low_N, cancer = high_Y + low_Y) %>% 
  select(study, control, cancer)

# Combine and aggregate the total n together
combine_stool <- adn_stool %>% full_join(crc_stool, by = "study") %>% 
  mutate(final_control = ifelse(is.na(control.x), invisible(control.y), invisible(control.x))) %>% 
  select(study, final_control, adenoma, cancer) %>% 
  rename(control = final_control)

# Read in the adenoma tissue files
adn_tissue <- read.csv("data/process/tables/alpha_adn_group_counts_tissue_summary.csv") %>% 
  filter(measure == "shannon") %>% 
  mutate(control = high_N + low_N, adenoma = high_Y + low_Y) %>% 
  select(study, control, adenoma)

# Read in the crc matched tissue files
crc_tissue_matched <- read_csv("data/process/tables/alpha_group_counts_matched_tissue_summary.csv") %>% 
  filter(measure == "shannon") %>% 
  mutate(control = high_N + low_N, cancer = high_Y + low_Y) %>% 
  select(study, control, cancer)

# Read in the crc unmatched tissue files
crc_tissue_unmatched <- read_csv("data/process/tables/alpha_group_counts_unmatched_tissue_summary.csv") %>% 
  filter(measure == "shannon") %>% 
  mutate(control = high_N + low_N, cancer = high_Y + low_Y) %>% 
  select(study, control, cancer)

# combine all the tissue samples together into a single file and aggregate
combine_tissue <- adn_tissue %>% full_join(crc_tissue_unmatched, by = "study") %>% 
  mutate(final_control = ifelse(is.na(control.x), invisible(control.y), invisible(control.x))) %>% 
  select(study, final_control, adenoma, cancer) %>% 
  rename(control = final_control) %>% 
  bind_rows(crc_tissue_matched) %>% 
  group_by(study) %>% 
  summarise(control = sum(control), adenoma = sum(adenoma), cancer = sum(cancer)) %>% 
  ungroup()


# Write out the tables
write_csv(combine_stool, "data/process/tables/stool_study_n_analyzed.csv")
write_csv(combine_tissue, "data/process/tables/tissue_study_n_analyzed.csv")



