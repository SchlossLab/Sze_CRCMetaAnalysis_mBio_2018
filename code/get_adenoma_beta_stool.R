### Get Stool Beta Diversity (Bray-Curtis)
### Look for overall group differences
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "lme4", "vegan"))

# Stool Only polyp sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("brim", "zeller", "baxter", "hale")



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
get_metadata_data <- function(i, sampleType){
  # i represents the data set
  # sampleType represents whether it is stool or tissue
  
  # Command that actually does the reading in and removes uneeded alpha metrics
  data_list <- read.csv(paste("data/process/tables/", i, "_", sampleType, "_",  
                              "transformed_data_alpha_raw_values.csv", sep = ""), 
                        header = T, stringsAsFactors = F) %>% 
    select(-sobs, -shannon, -shannoneven)
  
  # return to working environment the data list
  return(data_list)
}


# Function to make meta data the same order as that seen in the distance matrix
reorder_meta <- function(i, distanceList = distance_data, metaList = metadata){
  # i is the study
  # distanceList is defaulted to distance_data to make running mapply easier
  # metaList is defaulted to metadata to make running mapply easier
  
  # pulls the specific study distance matrix group order
  dist_names <- as.character(rownames(distanceList[[i]]))
  # aligns meta data by group order in distance matrix and creates a new is.cancer variable
  new_meta <- metaList[[i]] %>% slice(match(dist_names, as.character(group))) %>% 
    filter(disease != "cancer") %>% 
    mutate(is_adn =  factor(ifelse(disease == "polyp", 
                                      invisible("Y"), invisible("N")), levels = c("Y", "N")))
  # returns the modified metadata file
  return(new_meta)
}


# Function to align distance matrices with the metadata file
reorder_dist <- function(i, distanceList = distance_data, metaList = reordered_meta){
  # i is the study
  # distanceList is defaulted to distance_data to make running mapply easier
  # metaList is defaulted to reordered_meta so that the aligned metalist is used
  
  # grabs the samples in the meta file
  samples_to_keep <- as.character(metaList[[i]]$group)
  # subsets the distance matrix by the grabbed metafile samples
  new_distance_tbl <- distanceList[[i]][samples_to_keep, samples_to_keep]
  # returns the new updated distance matrix
  return(new_distance_tbl)
}


# Generate the PERMANOVA comparisons
make_adonis_test <- function(i, distanceList = reordered_dist, 
                             metaList = reordered_meta){
  # i is the study
  # distanceList is defaulted to reordered_dist so the subsetted dist file is used
  # metaList is defaulted to reordered_meta so the aligned metalist file is used
  
  # set the seed number so results don't keep changing
  set.seed(1234567)
  # run the permanova test using vegan for the selected data set
  temptest <- adonis(as.dist(distanceList[[i]]) ~ metaList[[i]]$is_adn, 
                     permutations = 9999)
  # pull only results of interest to be saved
  result_vector <- c(fstat = temptest$aov.tab$F.Model[1], r2 = temptest$aov.tab$R2[1], 
                     pvalue = temptest$aov.tab$`Pr(>F)`[1])
  # Print out progress
  print(paste("completed comparing ", i, sep = ""))
  # return the results of interest
  return(result_vector)
}


##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


# Read in distance data
distance_data <- mapply(get_distance, stool_sets, 
                        "data/process/", "braycurtis.0.03.lt.ave.dist")

# Read in data with needed metadata
metadata <- mapply(get_metadata_data, stool_sets, "stool", SIMPLIFY = F)

# Reorder the metadata to match distance row names
reordered_meta <- mapply(reorder_meta, stool_sets, SIMPLIFY = F)

# Reorder the distance matrix to match the the metadata length
reordered_dist <- mapply(reorder_dist, stool_sets, SIMPLIFY = F)

# Generate the total n in each group
total_compared <- lapply(reordered_meta, function(x) table(x$is_adn)) %>% 
  bind_cols() %>% t() %>% as.data.frame() %>% mutate(study = rownames(.)) %>% 
  rename(polyp_n = V1, control_n = V2)

# Get comparisons
beta_perm_results <- t(mapply(make_adonis_test, stool_sets)) %>% 
  as.data.frame() %>% mutate(study = rownames(.)) %>% 
  inner_join(total_compared, by = "study")


# Write out the data
write.csv(beta_perm_results, 
          "data/process/tables/beta_perm_adn_stool_summary.csv", row.names = F)




