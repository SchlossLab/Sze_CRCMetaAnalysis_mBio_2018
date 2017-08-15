### Get Stool Beta Diversity (Bray-Curtis)
### Look for overall group differences
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "lme4"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")


#### Need to do list
###### load data in
###### Use distance matrices to generate PERMANOVA values
###### Think of potential way to pool this information together
    ##### E.g. distance of centroids from each other...





