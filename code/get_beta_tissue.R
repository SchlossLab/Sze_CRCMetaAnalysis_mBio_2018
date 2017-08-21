### Get Tissue Beta Diversity (Bray-Curtis)
### Look for overall group differences matched vs unmatched
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "lme4", "vegan"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("dejea", "geng", "sana", "burns")

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


# Create function to read in distance data
get_distance <- function(i, path_to_file, file_ending){
  # i is the respective study
  # path_to_file is the path to the dir where file is stored
  # file_ending customizes the type of file that needs to be pulled
  
  # reads in the distance file based on variables provided
  tempDist <- read.dist(paste(path_to_file, i, "/", i, ".", file_ending, sep = ""))
  # returns the read in distance file
  return(tempDist)
}



# This function below is used as a proxy of the full metadata table since only 
# need a small component and it is contained in these files
get_metadata_data <- function(i,dataName){
  # i represents the data set
  # sampleType represents whether it is stool or tissue
  
  tempdata <- get(dataName)
  
  # Command that actually does the reading in and removes uneeded alpha metrics
  filtered_temp_data <- tempdata %>% filter(study == i) %>% 
    select(group, id, disease, study) %>% 
    mutate(id = as.character(id), 
           is_cancer = factor(ifelse(
             disease == "cancer", invisible("Y"), invisible("N")), levels = c("Y", "N")))
  
  # return to working environment the data list
  return(filtered_temp_data)
}


# Function to align distance matrices with the metadata file
reorder_dist <- function(i,  meta_name, distanceList = distance_data){
  # i is the study
  # distanceList is defaulted to distance_data to make running mapply easier
  # metaList is defaulted to reordered_meta so that the aligned metalist is used
  
  ref_data <- get(meta_name)
  
  # grabs the samples in the meta file
  samples_to_keep <- as.character(ref_data[[i]]$group)
  # subsets the distance matrix by the grabbed metafile samples
  new_distance_tbl <- distanceList[[i]][samples_to_keep, samples_to_keep]
  # returns the new updated distance matrix
  return(new_distance_tbl)
}


# Generate the PERMANOVA comparisons
make_adonis_test <- function(i, distanceList_name, 
                             metaList_name){
  # i is the study
  # distanceList is defaulted to reordered_dist so the subsetted dist file is used
  # metaList is defaulted to reordered_meta so the aligned metalist file is used
  
  distanceList <- get(distanceList_name)
  metaList <- get(metaList_name)
  
  # set the seed number so results don't keep changing
  set.seed(1234567)
  # run the permanova test using vegan for the selected data set
  temptest <- adonis(as.dist(distanceList[[i]]) ~ metaList[[i]]$is_cancer, 
                     permutations = 9999)
  # pull only results of interest to be saved
  result_vector <- c(fstat = temptest$aov.tab$F.Model[1], r2 = temptest$aov.tab$R2[1], 
                     pvalue = temptest$aov.tab$`Pr(>F)`[1])
  
  # return the results of interest
  return(result_vector)
}













# Read in distance data
distance_data <- mapply(get_distance, c(tissue_sets, both_sets), 
                        "data/process/", "braycurtis.0.03.lt.ave.dist")


# Read in data with needed metadata
unmatched_meta <- mapply(get_metadata_data, 
                         c("sana", "burns", both_sets), "tissue_unmatched", SIMPLIFY = F)


matched_meta <- mapply(get_metadata_data, 
                       c("dejea", "geng", "burns"), "tissue_matched", SIMPLIFY = F)


# Reorder the distance matrix to match the the metadata length
reordered_unmatched_dist <- mapply(reorder_dist, c("sana", "burns", both_sets), 
                                   "unmatched_meta", SIMPLIFY = F)

reordered_matched_dist <- mapply(reorder_dist, c("dejea", "geng", "burns"), 
                                 "matched_meta", SIMPLIFY = F)


# Get comparisons
beta_perm_unmatched_results <- t(mapply(make_adonis_test, c("sana", "burns", both_sets), 
                                        "reordered_unmatched_dist", "unmatched_meta")) %>% 
  as.data.frame() %>% mutate(study = rownames(.))


beta_perm_matched_results <- t(mapply(make_adonis_test, c("dejea", "geng", "burns"), 
                                        "reordered_matched_dist", "matched_meta")) %>% 
  as.data.frame() %>% mutate(study = rownames(.))




