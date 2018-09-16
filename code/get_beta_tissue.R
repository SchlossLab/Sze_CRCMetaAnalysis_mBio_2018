### Get Tissue Beta Diversity (Bray-Curtis)
### Look for overall group differences matched vs unmatched
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "lme4", "vegan"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu and Sana since it only has polyps and no cancer cases
tissue_sets <- c("dejea", "geng", "burns")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")


# Read in specific data tables to be used

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0)) %>% 
  filter(id != "3776") # remove polyp sample

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
  # dataName is a character string of the data file 
  
  # this grabs the data file from the global environment
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
  # meta_name is a character string of the meta file list of interest
  # distanceList is defaulted to distance_data to make running mapply easier
  
  # This grabs the meta file list from the global environment
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
  # distanceList_name is a chr string of the distance file list
  # metaList_name is a chr string of the meta file list
  
  # Grabs the distance and meta file list from the global environment
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



# Get specific bray-curtis distance between matched samples
get_bray_dist <- function(i, distanceList_name, metaList_name, 
                          controls, cases, getCont = FALSE, getCase = FALSE){
  # i is the study
  # distanceList_name is a chr string of the distance file list
  # metaList_name is a chr string of the meta file list
  # controls is what the control group is called
  # cases is what the case group is called
  # getCont controls whether we are collecting distances for controls only
  # getCase controls whether we are collecting distances for cases only
  
  # Sets up empty vectors to store distances
  c1_ids <- c() # control
  c2_ids <- c() # control
  can1_ids <- c() # case
  can2_ids <- c() # case
  # grabs the meta data file list from global environment
  metaList <- get(metaList_name)
  # grabs the sample ids for matched control and cases within the same person
  control_ids <- (metaList[[i]] %>% filter(disease == controls))[, "group"]
  case_ids <- (metaList[[i]] %>% filter(disease == cases))[, "group"]
  # if we want no control or case only distances go here otherwise...
  if(getCont == FALSE){
    # gets the bray distance for each specific pair of interest
    temp_values <- as.numeric(mapply(get_bray_value, i, case_ids, control_ids, distanceList_name, 
                                     USE.NAMES = F))
    # this is where we go if we want only control or case distances
  } else{
    # created a loop to add ids based on length of original vector of matched samples
    for(j in 1:length(control_ids)){
      # iteratively adds the ids with one less each id down it moves
      c1_ids <- c(c1_ids, rep(control_ids[j], (length(control_ids) - j)))
      can1_ids <- c(can1_ids, rep(case_ids[j], (length(case_ids) - j)))
      # do not want to add the last id since by then there are no unique comparisons
      if(j != length(control_ids)){
        # iteratively adds the ids with one less each id down it moves 
        c2_ids <- c(c2_ids, rev(control_ids)[1:(length(control_ids) - j)])
        can2_ids <- c(can2_ids, rev(case_ids)[1:(length(case_ids) - j)])
        # goes here when we hit the final length value
      } else {
        # simple saves the existing vector as itself
        c2_ids <- c2_ids
        can2_ids <- can2_ids
      }

    }
    # controls which distance values to grab: cases only or controls only
    if(getCase == FALSE){
      # gets the bray distance for each specific pair of interest for controls only
      temp_values <- as.numeric(mapply(get_bray_value, i, c1_ids, c2_ids, distanceList_name, 
                                       USE.NAMES = F))
    } else{
      # gets the bray distance for each specific pair of interest for cases only
      temp_values <- as.numeric(mapply(get_bray_value, i, can1_ids, can2_ids, distanceList_name, 
                                       USE.NAMES = F))
    }
    
  }
  # return the requested distances
  return(temp_values)
}


# Function to grab the needed data
get_bray_value <- function(i, row_value, col_value, 
                           distanceList_name){
  # i is the study
  # row_value is the first sample ID
  # col_value is the second sample ID
  # distanceList_name is a chr string of the distance file list
  
  # Grabs the distance file list of interest from the global environment
  distanceList <- get(distanceList_name)
  # grabs the specific distance value
  dist_value <- distanceList[[i]][row_value, col_value]
  # returns the specific distance value
  return(dist_value)
}


# Function to get wilcoxson p-value
run_wilcox <- function(i, cases, controls){
  # i is for study
  # cases are the distance value list for the first group
  # controls are the distance value list for the second group
  
  # this grabs the distance lists from the global environment
  cases <- get(cases)
  controls <- get(controls)
  
  # Get the respective values to be outputted
  tempData <- c(
  median_w_person = median(cases[[i]]), 
  q25_w_person = quantile(cases[[i]], probs = 0.25), 
  q75_w_person = quantile(cases[[i]], probs = 0.75), 
  median_bw_group = median(controls[[i]]), 
  q25_bw_group = quantile(controls[[i]], probs = 0.25), 
  q75_bw_group = quantile(controls[[i]], probs = 0.75), 
  pvalue = wilcox.test(cases[[i]], controls[[i]])$p.value)
  
  # returns the pvalue of the wilcoxson comparison between the two groups
  return(tempData)
  
}


##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


# Read in distance data
distance_data <- mapply(get_distance, c(tissue_sets, both_sets), 
                        "data/process/", "braycurtis.0.03.lt.ave.dist")


# Read in data with needed metadata
unmatched_meta <- mapply(get_metadata_data, 
                         c("burns", both_sets), "tissue_unmatched", SIMPLIFY = F)


matched_meta <- mapply(get_metadata_data, 
                       c("dejea", "geng", "burns"), "tissue_matched", SIMPLIFY = F)


# Reorder the distance matrix to match the the metadata length
reordered_unmatched_dist <- mapply(reorder_dist, c("burns", both_sets), 
                                   "unmatched_meta", SIMPLIFY = F)

reordered_matched_dist <- mapply(reorder_dist, c("dejea", "geng", "burns"), 
                                 "matched_meta", SIMPLIFY = F)

# Generate the total n in each group
unmatched_total_compared <- lapply(unmatched_meta, function(x) table(x$is_cancer)) %>% 
  bind_cols() %>% t() %>% as.data.frame() %>% mutate(study = rownames(.)) %>% 
  rename(crc_n = V1, not_crc_n = V2)

matched_total_compared <- lapply(matched_meta, function(x) table(x$is_cancer)) %>% 
  bind_cols() %>% t() %>% as.data.frame() %>% mutate(study = rownames(.)) %>% 
  rename(crc_n = V1, not_crc_n = V2)

# Get comparisons
beta_perm_unmatched_results <- t(mapply(make_adonis_test, c("burns", both_sets), 
                                        "reordered_unmatched_dist", "unmatched_meta")) %>% 
  as.data.frame() %>% mutate(study = rownames(.)) %>% 
  inner_join(unmatched_total_compared, by = "study")


beta_perm_matched_results <- t(mapply(make_adonis_test, c("dejea", "geng", "burns"), 
                                        "reordered_matched_dist", "matched_meta")) %>% 
  as.data.frame() %>% mutate(study = rownames(.)) %>% 
  inner_join(matched_total_compared, by = "study")

# Get vectors of matched and unmatched values

matched_bray_casetocontrol <- mapply(get_bray_dist, c("dejea", "geng", "burns"), 
                             "reordered_matched_dist", "matched_meta", "control", "cancer")

matched_bray_controls <- mapply(get_bray_dist, c("dejea", "geng", "burns"), 
                                "reordered_matched_dist", "matched_meta", "control", "cancer", 
                                getCont = TRUE, getCase = FALSE)

matched_bray_cases <- mapply(get_bray_dist, c("dejea", "geng", "burns"), 
                             "reordered_matched_dist", "matched_meta", "control", "cancer", 
                             getCont = TRUE, getCase = TRUE)


# Generate wilcoxson p-values
# this is essentially testing if the matched samples are closer to each other
# than different samples in the same disease group
bray_distance_matched_test_cont <- t(as.data.frame(mapply(run_wilcox, c("dejea", "geng", "burns"), 
               "matched_bray_casetocontrol", "matched_bray_controls", 
               SIMPLIFY = F))) %>% as.data.frame() %>% 
  mutate(study = rownames(.), bh = p.adjust(pvalue, method = "BH"))


bray_distance_matched_test_cases <- t(as.data.frame(mapply(run_wilcox, c("dejea", "geng", "burns"), 
                                          "matched_bray_casetocontrol", "matched_bray_cases", 
                                          SIMPLIFY = F))) %>% as.data.frame() %>% 
  mutate(study = rownames(.), bh = p.adjust(pvalue, method = "BH"))


# Write out the relevant data tables 
write.csv(beta_perm_unmatched_results, "data/process/tables/beta_perm_unmatched_tissue_summary.csv", 
          row.names = F)

write.csv(beta_perm_matched_results, "data/process/tables/beta_perm_matched_tissue_summary.csv", 
          row.names = F)

write.csv(bray_distance_matched_test_cont, 
          "data/process/tables/bray_matched_v_cont_tissue_summary.csv", 
          row.names = F)

write.csv(bray_distance_matched_test_cases, 
          "data/process/tables/bray_matched_v_cases_tissue_summary.csv", 
          row.names = F)



