### Get Stool Beta Diversity (Bray-Curtis)
### Look for overall group differences
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "lme4", "vegan"))

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen since there is only one case
both_sets <- c("flemer")


# Create function to read in distance data
get_distance <- function(i, path_to_file, file_ending){
  
  tempDist <- read.dist(paste(path_to_file, i, "/", i, ".", file_ending, sep = ""))
  
  return(tempDist)
}


# This function below is used as a proxy of the full metadata table since only 
# need a small component and it is contained in these files
get_metadata_data <- function(i, sampleType){
  # i represents the data set
  # sampleType represents whether it is stool or tissue
  
  # Command that actually does the reading in
  data_list <- read.csv(paste("data/process/tables/", i, "_", sampleType, "_",  
                              "transformed_data_alpha_raw_values.csv", sep = ""), 
                        header = T, stringsAsFactors = F) %>% 
    select(-sobs, -shannon, -shannoneven)
  
  # return to working environment the data list
  return(data_list)
}


# Function to make meta data the same order as that seen in the distance matrix
reorder_meta <- function(i, distanceList = distance_data, metaList = metadata){
  
  dist_names <- as.character(rownames(distanceList[[i]]))
  
  new_meta <- metaList[[i]] %>% slice(match(dist_names, as.character(group))) %>% 
    mutate(is_cancer =  factor(ifelse(disease == "cancer", 
                                      invisible("Y"), invisible("N")), levels = c("Y", "N")))

  return(new_meta)
}



# Function to align distance matrices
reorder_dist <- function(i, distanceList = distance_data, metaList = reordered_meta){
  
  samples_to_keep <- as.character(metaList[[i]]$group)
  
  new_distance_tbl <- distanceList[[i]][samples_to_keep, samples_to_keep]
  
  return(new_distance_tbl)
}


# Generate the PERMANOVA comparisons
make_adonis_test <- function(i, distanceList = reordered_dist, 
                             metaList = reordered_meta){
  
  set.seed(1234567)
  temptest <- adonis(as.dist(distanceList[[i]]) ~ metaList[[i]]$is_cancer, 
                     permutations = 9999)
  
  result_vector <- c(fstat = temptest$aov.tab$F.Model[1], r2 = temptest$aov.tab$R2[1], 
                     pvalue = temptest$aov.tab$`Pr(>F)`[1])
  
  
  return(result_vector)
  
}





# Read in distance data
distance_data <- mapply(get_distance, c(stool_sets, both_sets), 
                        "data/process/", "braycurtis.0.03.lt.ave.dist")

# Read in data with needed metadata
metadata <- mapply(get_metadata_data, c(stool_sets, both_sets), "stool", SIMPLIFY = F)

# Reorder the metadata to match distance row names
reordered_meta <- mapply(reorder_meta, c(stool_sets, both_sets), SIMPLIFY = F)

# Reorder the distance matrix to match the the metadata length
reordered_dist <- mapply(reorder_dist, c(stool_sets, both_sets), SIMPLIFY = F)


# Get comparisons
beta_perm_results <- t(mapply(make_adonis_test, c(stool_sets, both_sets))) %>% 
  as.data.frame() %>% mutate(study = rownames(.))


#### Need to do list
###### Think of potential way to pool this information together
    ##### E.g. distance of centroids from each other...





