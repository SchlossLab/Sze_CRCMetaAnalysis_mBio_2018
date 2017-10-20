### Code to measure AUCs of genus between studies within disease (Adenoma or Carcinoma)
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))



# Load needed data tables (adenoma -- stool)
adn_stool_studies <- c("baxter", "brim", "hale", "zeller")

adn_all_stool <- make_table(adn_stool_studies, "data/process/tables/", 
                            "adn_genus_stool_RF_full_", "adn_genus_stool_RF_select_", 
                            "_pvalue_summary.csv")

# Load needed data (adenoma -- tissue)
adn_tissue_studies <- c("flemer", "lu")

adn_all_tissue <- make_table(adn_tissue_studies, "data/process/tables/", 
                            "adn_genus_unmatched_tissue_RF_full_", 
                            "adn_genus_unmatched_tissue_RF_select_", 
                            "_pvalue_summary.csv")


# Load needed data (carcinoma -- stool)
crc_stool_studies <- c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller")

crc_all_stool <- make_table(crc_stool_studies, "data/process/tables/", 
                             "genus_stool_RF_full_", 
                             "genus_stool_RF_select_", 
                             "_pvalue_summary.csv")

# Load needed data (carcinoma -- tissue)
crc_matched_tissue_studies <- c("burns", "dejea", "geng")

crc_all_matched_tissue <- make_table(crc_matched_tissue_studies, "data/process/tables/", 
                            "genus_matched_tissue_RF_", 
                            "genus_matched_tissue_RF_select_", 
                            "_pvalue_summary.csv")

##############################################################################################
############################## List of functions to be used  #################################
##############################################################################################

# Function to generate the tables to be used
make_table <- function(studies, path_to_file, first_data_part_name, 
                       second_data_part_name, ending){
  
  tempFull <- sapply(studies, 
                     function(x) read_csv(paste(path_to_file, first_data_part_name, x, 
                                                ending, sep = "")) %>% 
                       mutate(model = x, model_type = "full"), simplify = F) %>% bind_rows()
  
  tempSelect <- sapply(studies, 
                       function(x) read_csv(paste(path_to_file, second_data_part_name, x, 
                                                  ending, sep = "")) %>% 
                         mutate(model = x, model_type = "select"), simplify = F) %>% bind_rows()
  
  tempALL <- tempFull %>% bind_rows(tempSelect)
  
  return(tempALL)
  
}


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