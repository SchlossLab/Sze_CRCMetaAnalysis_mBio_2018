### Code for the stool top 10 most important genera and OTUs Tax ID
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


crc_genera_occurances <- read_csv("data/process/tables/crc_RF_genera_stool_top10.csv")
adn_genera_occurances <- read_csv("data/process/tables/adn_RF_genera_stool_top10.csv")
crc_otu_occurances <- read_csv("data/process/tables/crc_RF_otu_stool_top10.csv")
adn_otu_occurances <- read_csv("data/process/tables/adn_RF_otu_stool_top10.csv")




##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################