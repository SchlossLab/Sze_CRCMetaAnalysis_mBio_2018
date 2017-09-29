### Pull, transform, and normalize 4 crc genera - adenoma
### Specifically analyze the genera tied to crc from previous research (tissue)
### One for unmatched and one for matched if possible
### Marc Sze

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "epiR", "metafor"))


# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("lu")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer")

# CRC genera of interest
crc_genera <- c("Fusobacterium", "Peptostreptococcus", "Porphyromonas", "Parvimonas")



##############################################################################################
############### List the functions to be used for analysis ###################################
##############################################################################################















##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################







