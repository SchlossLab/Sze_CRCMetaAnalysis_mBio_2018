### Generate total n used for each study
### Use the RR based count to generate this data
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse"))

# Read in needed data tables

adn_stool <- read_csv("data/process/tables/alpha_adn_group_counts_summary.csv") %>% 
  filter(measure == "sobs") %>% 
  mutate(control = high_N + low_N)
  
  
crc_stool <- read_csv("data/process/tables/alpha_group_counts_summary.csv") %>% 
  filter(measure == "shannon") %>% 
  mutate(cancer = high_Y + low_Y)

adn_tissue <- read.csv("data/process/tables/alpha_adn_group_counts_tissue_summary.csv") %>% 
  filter(measure == "shannon")

crc_tissue_matched <- read_csv("data/process/tables/alpha_group_counts_matched_tissue_summary.csv") %>% 
  filter(measure == "shannon")

crc_tissue_unmatched <- read_csv("data/process/tables/alpha_group_counts_tissue_summary.csv") %>% 
  filter(measure == "shannon")



