### Get the dmm groups - tissue
### Merge the group data with relevant metadata
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "vegan", "foreach", "doParallel"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "flemer")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")