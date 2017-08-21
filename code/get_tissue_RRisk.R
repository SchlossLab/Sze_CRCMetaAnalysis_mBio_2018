### Get Tissue relative risk comparisons
### One for unmatched and one for matched if possible
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "epiR", "metafor"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
# Remove dejea since there is only one without cancer
tissue_sets <- c("sana", "burns", "geng")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")


# Read in specific data tables to be used

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0))

tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  mutate(disease = gsub("adenoma", "polyp", disease))

# Remove polyp only group
no_p_tissue_unmatched <- tissue_unmatched %>% filter(study != "lu" & study != "dejea")

# Function to run the analysis
analyze_study <- function(i, var_of_int, group_column, dataset = no_p_tissue_unmatched){
  
  working_data <- dataset %>% filter(study == i)
  thresholds <- apply(select(working_data, 
                             one_of("r_sobs", "r_shannon", "r_shannoneven")), 2, 
                      function(x) median(x))
  
  # Generate the disease group vector
  is_cancer <- factor(ifelse(working_data[, group_column] == "cancer", 
                             invisible("Y"), invisible("N")), levels = c("Y", "N"))
  
  # Runs the code to generate high/low calls for the alpha metrics used based on median
  highs_lows <- mapply(create_high_low, i, c("r_sobs", "r_shannon", "r_shannoneven"), thresholds, 
                       "disease", SIMPLIFY = F)
  names(highs_lows) <- c("sobs", "shannon", "shannoneven") # forces names for the list
  
  return(is_cancer)
  
}


# Function that creates the needed high/low columns
create_high_low <- function(i, var_of_interest, threshold, grouping, 
                            dataset = tissue_unmatched){
  # i is the study
  # var_of_interest is the alpha metrics being used
  # threshold is the vector of median values for alpha measures of interest
  # grouping is the name of the case/control column
  # dataset is default to the tissue_unmatched to allow for mapply to work
  
  working_data <- dataset %>% filter(study == i)
  
  # create a vector with high/low versus the median value provided
  high_low <- factor(ifelse(working_data[, var_of_interest] <= threshold, 
                            invisible("low"), invisible("high")), levels = c("low", "high"))
  # Returns the vector of high/low calls
  return(high_low)
}










