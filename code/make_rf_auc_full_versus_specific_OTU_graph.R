### Code to measure AUCs of full data set and select OTU dataset
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Load needed data tables (adenoma)
adn_tissue <- read_csv("data/process/tables/adn_tissue_rf_otu_random_comparison_summary.csv") %>% 
  mutate(type = c("matched", "unmatched"))

adn_all_stool <- read_csv("data/process/tables/adn_stool_rf_otu_random_comparison_summary.csv") %>% 

# Load in needed data tables (carcinoma)



##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

### Study Colors by Viridis 
##  library(scales) 
## show_col(viridis_pal()(16))
# flemer - #440154FF
# lu - #FDE725FF
# burns - #453581FF
# chen - #3D4D8AFF
# sana - #1F998AFF
# dejea - #2B748EFF
# geng - #CBE11EFF
# brim - #34618DFF
# zeller - #FDE725FF
# baxter - #481D6FFF
# hale - #67CC5CFF
# wang - #97D83FFF
# weir - #24878EFF
# ahn - #40BC72FF





##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

