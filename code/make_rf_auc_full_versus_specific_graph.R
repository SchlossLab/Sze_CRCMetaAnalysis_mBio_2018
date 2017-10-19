### Code to measure AUCs of full data set and select genera dataset
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Load needed data tables (adenoma)
adn_tissue_matched <- 
  read_csv("data/process/tables/adn_genus_matched_tissue_RF_fullvsselect_pvalue_summary.csv")
adn_tissue_unmatched <- 
  read_csv("data/process/tables/adn_genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv")

adn_all_stool <- read_csv("data/process/tables/adn_genus_stool_RF_fullvsselect_pvalue_summary.csv")

# Load in needed data tables (carcinoma)
crc_tissue_matched <- 
  read_csv("data/process/tables/genus_matched_tissue_RF_fullvsselect_pvalue_summary.csv")
crc_tissue_unmatched <- 
  read_csv("data/process/tables/genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv")

crc_all_stool <- read_csv("data/process/tables/genus_stool_RF_fullvsselect_pvalue_summary.csv")
