### Code to create Specific Genus RR pooled graphs
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))

# Load needed data tables (adenoma)
adn_all_stool <- read_csv("data/process/tables/adn_select_genus_RR_stool_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_RR_stool_ind_results.csv"))

adn_all_tissue <- read_csv("data/process/tables/adn_select_genus_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_RR_tissue_ind_results.csv"))

# Load in needed data tables (carcinoma)
crc_all_stool <- read_csv("data/process/tables/select_genus_RR_stool_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/select_genus_RR_stool_ind_results.csv"))


crc_all_tissue <- read_csv("data/process/tables/select_genus_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/select_genus_RR_tissue_ind_results.csv"))



