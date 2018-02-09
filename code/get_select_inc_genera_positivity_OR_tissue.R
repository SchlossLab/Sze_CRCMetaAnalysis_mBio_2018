### Test each increasing positive amount against 0 for RR of CRC or adenoma
### Tissue
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "epiR", "metafor"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("dejea", "geng", "sana", "burns")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")

# All tissue studies
all_studies <- c("burns", "chen", "flemer", "sana", "dejea", "geng")

# CRC genera of interest
crc_genera <- c("Fusobacterium", "Peptostreptococcus", "Porphyromonas", "Parvimonas")


##############################################################################################
############### List of Programs to make the analysis Run ####################################
##############################################################################################

# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_data <- function(i, metadata_table){
  # i is the study of interest
  
  # gets original sample names
  sample_names <- (get_file(i, "data/process/", "_genera_shared.csv") %>% 
                     mutate(sample_names = rownames(.)) %>% select(sample_names))[, "sample_names"]
  # grabs subsampled data and assigns rownames from sample names to table
  sub_genera_data <- get_file(i, "data/process/", "_subsample_genera.csv", rows_present = F, 
                              vec_of_rownames = sample_names) %>% 
    as.data.frame() %>% mutate(sample_ID = rownames(.)) %>% 
    select(sample_ID, everything())
  # grabs the meta data and transforms polyp to control (polyp/control vs cancer) 
  study_meta <- metadata_table %>% filter(study == i) %>% rename(sampleID = group)
  
  # conditional that checks for whether length of rows of meta data is smaller
  if(length(rownames(study_meta)) < length(rownames(sub_genera_data))){
    # grab only the samples in the meta data file for down stream analysis
    sub_genera_data <- sub_genera_data %>% slice(match(study_meta$sampleID, sample_ID))
    
  } else{
    # grab only files in the data file for analysis
    study_meta <- study_meta %>% slice(match(sub_genera_data$sample_ID, sampleID))
  }
  # re assigns the rown names while removing the extra column used for sorting
  sample_names <- sub_genera_data$sample_ID
  # creates a list file with both data sets
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta)
  # returns the combined list file
  return(dataList)
  
}


# Function to grab only the genera file and pull specific genus from it
get_specific_genera <- function(i, genera_to_get, table_name, meta_name, 
                                dataList){
  # i represents the study
  # genera_to_get represents the genera of interest to pull specifically
  # table_name represents the table the has all the genera data (subsampled)
  # meta_name represents the metadata name that stores the relevent meta data
  # dataList is a list with both meta data and sub sampled genus data
  
  # grab the specific genera and merge with the meta data file
  tempData <- dataList[[i]][[table_name]] %>% 
    select(sample_ID, one_of(genera_to_get)) %>% 
    rename(sampleID = sample_ID) %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    inner_join((dataList[[i]][[meta_name]] %>% 
                  mutate(sampleID = as.character(sampleID), 
                         disease = ifelse(disease == "normal", 
                                          invisible("control"), invisible(disease)))), 
               by = "sampleID") %>% 
    mutate(disease2 = ifelse(disease != "cancer", invisible("control"), invisible(disease))) %>% 
    as.data.frame()
  
  # Print output when merged completed
  print(paste("Finished pulling and merging files for", i, "data sets"))
  # Return the newly transformed data
  return(tempData)
  
}


# Function to generate total positives and overalls
get_select_group_totals <- function(i, select_genera, dataList){
  
  tempData <- dataList[[i]] %>% 
    mutate_at(select_genera, 
              function(x) ifelse(x > 0, invisible(1), invisible(0))) %>% 
    mutate(all_four = rowSums(.[, select_genera]))
  
  return(dataList[[i]] %>% mutate(all_four = tempData$all_four))
  
}

# Function to generate the counts
generate_counts <- function(i, total_genera_present, dataList){
  
  tempData <- dataList[[i]]
  temp_two_by_two <- c(high_N = length(rownames(filter(tempData, 
                                  disease == "control" & all_four == total_genera_present))), 
                       high_Y = length(rownames(filter(tempData, 
                                  disease == "cancer" & all_four == total_genera_present))), 
                       low_N = length(rownames(filter(tempData, 
                                  disease == "control" & all_four == 0))), 
                       low_Y = length(rownames(filter(tempData, 
                                   disease == "cancer" & all_four == 0)))) 
  
  
  return(temp_two_by_two)
}


# Function to run test for selected alpha measure
run_pooled <- function(i, name_vector, dataList){
  # alpha_d is the alpha measure of interest
  # dataset is defaulted to ind_counts_data
  
  # select only the relevent alpha measures
  test_data <- dataList[[i]]
  
  # Run the actual pooled test
  rr_pooled_test <- rma(ai = high_Y, bi = high_N, 
                        ci = low_Y, di = low_N, data = test_data, 
                        measure = "RR", method = "REML")
  # Store a vector of the important results of interest
  results <- c(exp(c(rr = rr_pooled_test$b[[1, 1]], ci_lb = rr_pooled_test$ci.lb, 
                     ci_ub=rr_pooled_test$ci.ub)), pvalue = rr_pooled_test$pval, 
               measure = name_vector[i])
  # returns the vector of results
  return(results)
  
}






##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# Read in specific data tables to be used

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0)) %>% 
  select(-sobs, -shannon, -shannoneven, -r_sobs, -r_shannon, -r_shannoneven)

tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  mutate(disease = gsub("adenoma", "polyp", disease)) %>% 
  select(-sobs, -shannon, -shannoneven, -r_sobs, -r_shannon, -r_shannoneven)


# Manually remove the one matched polyp sample in dejea study
tissue_matched <- tissue_matched %>% filter(id != 3776)

# Create new column where polyp is assigned to control for tissue unmatched group
tissue_unmatched <- tissue_unmatched %>% 
  mutate(disease2 = ifelse(disease == "polyp", invisible("control"), invisible(disease)))


# Remove polyp only group
no_p_tissue_matched <- tissue_matched %>% filter(study != "lu")
no_p_tissue_unmatched <- tissue_unmatched %>% filter(study != "lu")
no_p_tissue <- tissue_matched %>% bind_rows(tissue_unmatched) %>% filter(study != "lu")

# Generate RR and data tables for every study
ind_data <- sapply(all_studies,  
                   function(x) get_data(x, no_p_tissue), simplify = F)


specific_genera_list <- sapply(all_studies, 
                               function(x) get_specific_genera(x, crc_genera, 
                                                               "sub_genera_data", "study_meta", 
                                                               ind_data), simplify = F)
# Get specific grouping with all the big 4 considered
mod_specific_genera_list <- sapply(all_studies, 
                                   function(x) 
                                     get_select_group_totals(x, crc_genera, specific_genera_list), 
                                   simplify = F)

# Generate the counts from 1 to 4 for the RR
overall_counts <- list()
select_counts <- c(one_counts = 1, two_counts = 2, three_counts = 3, four_counts = 4) 
count_names <- list("one_counts", "two_counts", "three_counts", "four_counts")


test <- sapply(count_names, 
               function(x) 
                 overall_counts[[x]] = as.data.frame(t(sapply(all_studies, 
                    function(y) 
                      generate_counts(y, select_counts[x], mod_specific_genera_list)))), simplify = F)


pooled_results <- sapply(c(1:4), 
                         function(x) 
                           run_pooled(x, count_names, test), simplify = F) %>% bind_rows()


write_csv(pooled_results, "data/process/tables/select_genus_inc_4_tissue.csv")


