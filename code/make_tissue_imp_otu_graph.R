### Code for the tissue top 10 most important genera and OTUs Tax ID
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis", "stringr"))


# Read in Needed data
write_csv(crc_unmatched_genera_occurances, "data/process/tables/crc_RF_genera_unmatched_tissue_top10.csv")
write_csv(crc_matched_genera_occurances, "data/process/tables/adn_RF_genera_matched_tissue_top10.csv")
write_csv(adn_genera_occurances, "data/process/tables/adn_RF_genera_tissue_top10.csv")
write_csv(crc_unmatched_otu_occurances, "data/process/tables/crc_RF_otu_unmatched_tissue_top10.csv")
write_csv(crc_matched_otu_occurances, "data/process/tables/adn_RF_otu_matched_tissue_top10.csv")
write_csv(adn_otu_occurances, "data/process/tables/adn_RF_otu_tissue_top10.csv")


adn_tissue_sets <- 2
crc_unmatched_sets <- 4
crc_matched_sets <- 3

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################





















##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


stool_graph <- grid.arrange(adn_genera, crc_genera, adn_otu, crc_otu)


ggsave("results/figures/Figure6.pdf", 
       stool_graph, width = 11, height = 8, dpi = 300)