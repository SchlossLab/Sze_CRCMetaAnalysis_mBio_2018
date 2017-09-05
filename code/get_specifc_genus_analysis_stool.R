### Pull, transform, and normalize 4 crc genera 
### Specifically analyze the genera tied to crc from previous research
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "epiR", "metafor"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")

# CRC genera of interest
crc_genera <- c("Fusobacterium", "Peptostreptococcus", "Porphyromonas", "Parvimonas")



# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_data <- function(i){
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
  study_meta <- get_file(i, "data/process/", ".metadata", rows_present = F,  
                         "stool", metadata = T)
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
                                dataList = stool_study_data){
  # i represents the study
  # genera_to_get represents the genera of interest to pull specifically
  # table_name represents the table the has all the genera data (subsampled)
  # meta_name represents the metadata name that stores the relevent meta data
  # dataList is defaulted to stool_study_data for convience
  
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

########################################################################################
### not possible to power transform to a normal distribution. Too much 0 weighting #####
### Main other possibility is to use RR and above/below median value               #####
########################################################################################

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
  highs_lows <- mapply(create_high_low, i, thresholds, 
                       vec_of_int, group_column, SIMPLIFY = F)
  names(highs_lows) <- vec_of_int # forces names for the list
  # Obtains the individual relative risk and CI for each study
  obtained_rr <- lapply(highs_lows, 
                        function(x) run_rr(high_low_vector = x, disease_vector = is_cancer)) 
  return(obtained_rr) # returns a list that stores all the table counts and RR data
}


# Function that creates the needed high/low columns
create_high_low <- function(i, threshold, var_of_interest, grouping, 
                            dataset = specific_genera_list){
  # i is the study
  # threshold is the vector of median values for alpha measures of interest
  # var_of_interest is the genera being used
  # grouping is the name of the case/control column
  # dataset is default to the stool_data list to allow for mapply to work
  
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
pull_data <- function(var_of_int, i, result, datalist =  test_ind_RR){
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
  pulled_data <- mapply(pull_data, vec_of_interest, 
                        i, result, SIMPLIFY = F) %>% bind_rows()
  # returns a nice data table
  return(pulled_data)
}


# Function to run test for selected alpha measure
run_pooled <- function(alpha_d, dataset = ind_counts_data){
  # alpha_d is the alpha measure of interest
  # dataset is defaulted to ind_counts_data
  
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

# reads in all the stool data into one list
stool_study_data <- mapply(get_data, c(stool_sets, "flemer"), SIMPLIFY = F)

# pull the specific genera of interest and merge with the meta data
specific_genera_list <- sapply(c(stool_sets, "flemer"), 
               function(x) get_specific_genera(x, crc_genera, 
                                               "sub_genera_data", "study_meta"))

# Generate the RR for each respective study for each genus of interest
# Return both counts and results
test_ind_RR <- sapply(c(stool_sets, "flemer"), 
               function(x) analyze_study(x, "disease2",crc_genera, specific_genera_list), 
               simplify = F)

# Store the results from the individual testing here
ind_RR_data <- sapply(c(stool_sets, "flemer"), 
                        function(x) make_list(x, crc_genera, "test_values", test_ind_RR), 
                      simplify = F) %>% bind_rows()

# Store the counts and rearrange the table to be used in the pooled analysis
ind_counts_data <- sapply(c(stool_sets, "flemer"), 
               function(x) make_list(x, crc_genera, "data_tbl", test_ind_RR), simplify = F) %>% 
  bind_rows() %>% unite(group, high_low_vector, disease_vector, sep = "_") %>% 
  spread(group, Freq)

# Run the pooled analysis for each respective genera of interest
pooled_results <- t(mapply(run_pooled, crc_genera, USE.NAMES = F)) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  mutate_at(c("rr", "ci_lb", "ci_ub", "pvalue"), as.numeric)






