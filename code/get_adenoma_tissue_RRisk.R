### Get Tissue relative risk comparisons - adenoma
### One for unmatched and one for matched if possible
### Marc Sze

# For this analysis combined both unmatched and matched together to boost n and power
# Can label which studies had majority matched to see if that made a difference at all.

# Load in needed functions and libraries
source('code/functions.R')

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












