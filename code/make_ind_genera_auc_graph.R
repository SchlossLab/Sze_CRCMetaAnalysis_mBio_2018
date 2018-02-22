### Code to create alpha RR graph Figure 2
### Stool specific measures only
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Load in needed data tables

stool_data <- read_csv("data/process/tables/ind_genera_auc_stool.csv")

combined_stool <- stool_data %>% group_by(taxa) %>% 
  summarise(auc = median(auc)) %>% 
  mutate(study = "median") %>% ungroup() %>% 
  bind_rows(stool_data)
  

unmatched_tissue_data <- read_csv("data/process/tables/ind_genera_auc_unmatched_tissue.csv")

combined_unmatched_tissue <- unmatched_tissue_data %>% group_by(taxa) %>% 
  summarise(auc = median(auc)) %>% 
  mutate(study = "median") %>% ungroup() %>% 
  bind_rows(unmatched_tissue_data)

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################
















##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

stool_alpha_RR <- grid.arrange(adn_stool_graph, crc_stool_graph)

ggsave("results/figures/Figure1.pdf", stool_alpha_RR, width = 8.5, height = 7, dpi = 300)