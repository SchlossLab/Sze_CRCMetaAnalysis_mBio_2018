### Pull, transform, and normalize 4 crc genera 
### Specifically analyze the genera tied to crc from previous research (tissue)
### One for unmatched and one for matched if possible
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
    mutate(all_four = rowSums(.[, select_genera]), 
           total_four = rowSums(dataList[[i]][, select_genera]))
  
  return(dataList[[i]] %>% mutate(all_four = tempData$all_four, total_four = tempData$total_four))
  
}



# Analyze the data with respect to high low column table
analyze_study <- function(i, group_column, vec_of_int, dataset){
  # i represents the study 
  # group_column represents what the case/control column is
  # vec_of_int is a vector with genera of interest
  # dataset is the list of combined genera of interest and meta data
  
  # Vector of whether sample was cancer or not
  is_cancer <- factor(ifelse(dataset[[i]][, group_column] == "cancer", 
                             invisible("Y"), invisible("N")), levels = c("Y", "N"))
  # Vector of median values 
  thresholds <- apply(select(dataset[[i]], one_of(vec_of_int)), 2, 
                      function(x) median(x))
  # Runs the code to generate high/low calls for the alpha metrics used based on median
  highs_lows <- sapply(vec_of_int, 
                       function(x) create_high_low(i = i, thresholds, 
                                                   x, group_column, 
                                                   dataset = dataset), simplify = F)

  names(highs_lows) <- vec_of_int # forces names for the list
  # Obtains the individual relative risk and CI for each study
  obtained_rr <- lapply(highs_lows, 
                        function(x) run_rr(high_low_vector = x, disease_vector = is_cancer)) 
  return(obtained_rr) # returns a list that stores all the table counts and RR data
}


# Function that creates the needed high/low columns
create_high_low <- function(i, threshold, var_of_interest, grouping, 
                            dataset){
  # i is the study
  # threshold is the vector of median values for alpha measures of interest
  # var_of_interest is the genera being used
  # grouping is the name of the case/control column
  # dataset is the data list you want to work with
  
  # get specific data table of interest based on study
  select_data <- dataset[[i]]
  
  # create a vector with high/low versus the median value provided
  high_low <- factor(ifelse(select_data[, var_of_interest] <= threshold, 
                            invisible("low"), invisible("high")), levels = c("high", "low"))
  # Returns the vector of high/low calls
  return(high_low)
}


# Function that runs relative risk test on single variable
run_rr <- function(high_low_vector, disease_vector){
  # high_low_vector is the respective call columns from high_low for a specific alpha measure
  # disease_vector is the "is_cancer" vector is case/control info
  
  # Creates a 2x2 table of counts
  contingency <- table(high_low_vector, disease_vector)
  
  # check if there are 0 counts
  check_values <- as.vector(contingency)
  
  # runs the RR test based on the obtained 2x2 table
  test <- epi.2by2(contingency, method="cohort.count")
  # Pull only specific information from the stored list in "test"
  test_values <- cbind(test$massoc$RR.strata.score, 
                       pvalue = test$massoc$chisq.strata$p.value)
  # store both the obtained raw counts and the resulting RR with pvalue
  combined_data <- list(data_tbl = contingency, test_values = test_values)
  
  # Returns a list with all information needed for downstream analysis
  return(combined_data)
}


# Function to seperate out the table data from the individual analysis data
pull_data <- function(var_of_int, i, result, datalist){
  # var_of_int is the alpha measures used e.g. "sobs"
  # i is the study
  # result is the type of data we want either "test_values" or "data_tbl"
  # datalist is defaulted to test_ind_RR to make it easier to work with mapply
  
  # Pull the needed data and add identifiers
  tempData <- datalist[[i]][[var_of_int]][[result]] %>% as.data.frame() %>% 
    mutate(measure = var_of_int, study = i)
  # return the pulled data
  return(tempData)
  
}


# A control function to direct final table creation for ind RR analysis
make_list <- function(i, vec_of_interest, result, datalist){
  # i is the study
  # result is they type of data being pulled "test_values" or "data_tbl"
  # vec_of_interest is the genera to be analyzed
  # datalist is defaulted to ind_study_data to make it easier to work with mapply
  
  # runs the function iteratively to collect the specific data
  pulled_data <- sapply(vec_of_interest, 
                        function(x) pull_data(i = i, x, 
                                              result = result, 
                                              datalist = datalist), simplify = F) %>% 
    bind_rows()
  
    #pulled_data <- mapply(pull_data, vec_of_interest, 
  #                      i, result, SIMPLIFY = F) %>% bind_rows()
  # returns a nice data table
  return(pulled_data)
}


# Function to run test for selected alpha measure
run_pooled <- function(alpha_d, dataset){
  # alpha_d is the alpha measure of interest
  # dataset is the individual count data 
  
  # select only the relevent alpha measures
  test_data <- dataset %>% filter(measure == alpha_d)
  
  # Run the actual pooled test
  rr_pooled_test <- rma(ai = high_Y, bi = high_N, 
                        ci = low_Y, di = low_N, data = test_data, 
                        measure = "RR", method = "REML")
  # Store a vector of the important results of interest
  results <- c(exp(c(rr = rr_pooled_test$b[[1, 1]], ci_lb = rr_pooled_test$ci.lb, 
                     ci_ub=rr_pooled_test$ci.ub)), pvalue = rr_pooled_test$pval, 
               measure = alpha_d)
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
ind_matched_data <- sapply(c("burns", "dejea", "geng"), 
                           function(x) get_data(x, no_p_tissue_matched), simplify = F)

ind_unmatched_data <- sapply(c("burns", "dejea", "sana", both_sets),  
                             function(x) get_data(x, no_p_tissue_unmatched), simplify = F)

ind_data <- sapply(all_studies,  
                   function(x) get_data(x, no_p_tissue), simplify = F)


# pull the specific genera of interest and merge with the meta data
matched_specific_genera_list <- sapply(
  c("burns", "dejea", "geng"), 
  function(x) get_specific_genera(x, crc_genera, 
                                  "sub_genera_data", "study_meta", 
                                  ind_matched_data), simplify = F)

unmatched_specific_genera_list <- sapply(
  c("burns", "sana", both_sets), 
  function(x) get_specific_genera(x, crc_genera, 
                                  "sub_genera_data", "study_meta", 
                                  ind_unmatched_data), simplify = F)


specific_genera_list <- sapply(all_studies, 
  function(x) get_specific_genera(x, crc_genera, 
                                  "sub_genera_data", "study_meta", 
                                  ind_data), simplify = F)


# Get specific grouping with all the big 4 considered
mod_matched_specific_genera_list <- sapply(c("burns", "dejea", "geng"), 
    function(x) get_select_group_totals(x, crc_genera, matched_specific_genera_list), simplify = F)

mod_unmatched_specific_genera_list <- sapply(c("burns", "sana", both_sets), 
    function(x) get_select_group_totals(x, crc_genera, unmatched_specific_genera_list), simplify = F)

mod_specific_genera_list <- sapply(all_studies, 
    function(x) get_select_group_totals(x, crc_genera, specific_genera_list), simplify = F)

# Generate the RR for each respective study for each genus of interest
# Return both counts and results
matched_test_ind_RR <- sapply(
  c("burns", "dejea", "geng"), 
  function(x) analyze_study(x, "disease", 
                            c(crc_genera, "all_four", "total_four"), 
                            mod_matched_specific_genera_list), simplify = F)


unmatched_test_ind_RR <- sapply(
  c("burns", "sana", both_sets), 
  function(x) analyze_study(x, "disease2", 
                            c(crc_genera, "all_four", "total_four"), 
                            mod_unmatched_specific_genera_list), simplify = F)

test_ind_RR <- sapply(all_studies, 
  function(x) analyze_study(x, "disease2", 
                            c(crc_genera, "all_four", "total_four"), 
                            mod_specific_genera_list), simplify = F)


# Store the results from the individual testing here
matched_RR_data <- sapply(c("burns", "dejea", "geng"), 
                      function(x) make_list(x, c(crc_genera, "all_four", "total_four"), 
                                            "test_values", matched_test_ind_RR), 
                      simplify = F) %>% bind_rows()

unmatched_RR_data <- sapply(c("burns", "sana", both_sets), 
                          function(x) make_list(x, c(crc_genera, "all_four", "total_four"), 
                                                "test_values", unmatched_test_ind_RR), 
                          simplify = F) %>% bind_rows()


RR_data <- sapply(all_studies, 
                            function(x) make_list(x, c(crc_genera, "all_four", "total_four"), 
                                                  "test_values", test_ind_RR), 
                            simplify = F) %>% bind_rows()

# Store the counts and rearrange the table to be used in the pooled analysis
matched_counts_data <- sapply(c("burns", "dejea", "geng"), 
                          function(x) make_list(x, c(crc_genera, "all_four", "total_four"), "data_tbl", 
                                                matched_test_ind_RR), simplify = F) %>% 
  bind_rows() %>% unite(group, high_low_vector, disease_vector, sep = "_") %>% 
  spread(group, Freq)



unmatched_counts_data <- sapply(c("burns", "sana", both_sets), 
                              function(x) make_list(x, c(crc_genera, "all_four", "total_four"), "data_tbl", 
                                                    unmatched_test_ind_RR), simplify = F) %>% 
  bind_rows() %>% unite(group, high_low_vector, disease_vector, sep = "_") %>% 
  spread(group, Freq)


counts_data <- sapply(all_studies, 
                                function(x) make_list(x, c(crc_genera, "all_four", "total_four"), 
                                                      "data_tbl", test_ind_RR), simplify = F) %>% 
  bind_rows() %>% unite(group, high_low_vector, disease_vector, sep = "_") %>% 
  spread(group, Freq)


# Run the pooled analysis for each respective genera of interest
matched_pooled_results <- t(sapply(c(crc_genera, "all_four", "total_four"), 
                         function(x) run_pooled(x, matched_counts_data))) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  mutate_at(c("rr", "ci_lb", "ci_ub", "pvalue"), as.numeric)


unmatched_pooled_results <- t(sapply(c(crc_genera, "all_four", "total_four"), 
                                   function(x) run_pooled(x, unmatched_counts_data))) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  mutate_at(c("rr", "ci_lb", "ci_ub", "pvalue"), as.numeric)


pooled_results <- t(sapply(c(crc_genera, "all_four", "total_four"), 
                                     function(x) run_pooled(x, counts_data))) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  mutate_at(c("rr", "ci_lb", "ci_ub", "pvalue"), as.numeric)

# Write out the important tables
write.csv(matched_counts_data, 
          "data/process/tables/select_genus_matched_tissue_group_counts_summary.csv", 
          row.names = F)
write.csv(matched_RR_data, "data/process/tables/select_genus_RR_matched_tissue_ind_results.csv", 
          row.names = F)
write.csv(matched_pooled_results, "data/process/tables/select_genus_RR_matched_tissue_composite.csv", 
          row.names = F)


write.csv(unmatched_counts_data, 
          "data/process/tables/select_genus_unmatched_tissue_group_counts_summary.csv", 
          row.names = F)
write.csv(unmatched_RR_data, "data/process/tables/select_genus_RR_unmatched_tissue_ind_results.csv", 
          row.names = F)
write.csv(unmatched_pooled_results, 
          "data/process/tables/select_genus_RR_unmatched_tissue_composite.csv", 
          row.names = F)

write.csv(counts_data, 
          "data/process/tables/select_genus_tissue_group_counts_summary.csv", 
          row.names = F)
write.csv(RR_data, "data/process/tables/select_genus_RR_tissue_ind_results.csv", 
          row.names = F)
write.csv(pooled_results, 
          "data/process/tables/select_genus_RR_tissue_composite.csv", 
          row.names = F)



