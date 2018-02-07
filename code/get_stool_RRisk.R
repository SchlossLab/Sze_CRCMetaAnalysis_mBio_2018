### Get Stool relative risk comparisons
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "epiR", "metafor"))


# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen since there is only one case
both_sets <- c("flemer")


# Function to read in the respective transformed data tables for each study
# These should have been previously power transformed to correct skew
get_data <- function(i, sampleType){
  # i represents the data set
  # sampleType represents whether it is stool or tissue
  
  # Command that actually does the reading in
  data_list <- read_csv(paste("data/process/tables/", i, "_", sampleType, "_",  
                              "alpha_raw_values.csv", sep = "")) %>% filter(disease != "polyp")
  
  # return to working environment the data list
  return(data_list)
}


# Analyze the data with respect to high low column table
analyze_study <- function(i, group_column, dataset = stool_data){
  # i represents the study 
  # group_column represents what the case/control column is
  # to work with mapply dataset is defaulted to read in list of data called stool_data
  
  # Vector of whether sample was cancer or not
  is_cancer <- factor(ifelse(dataset[[i]][, group_column] == "cancer", 
                      invisible("Y"), invisible("N")), levels = c("Y", "N"))
  # Vector of median values 
  thresholds <- apply(select(dataset[[i]], one_of("sobs", "shannon", "shannoneven")), 2, 
                      function(x) median(x))
  # Runs the code to generate high/low calls for the alpha metrics used based on median
  highs_lows <- mapply(create_high_low, i, thresholds, 
                       c("sobs", "shannon", "shannoneven"), "disease", SIMPLIFY = F)
  names(highs_lows) <- c("sobs", "shannon", "shannoneven") # forces names for the list
  # Obtains the individual relative risk and CI for each study
  obtained_rr <- lapply(highs_lows, 
                        function(x) run_rr(high_low_vector = x, disease_vector = is_cancer)) 
  return(obtained_rr) # returns a list that stores all the table counts and RR data
}


# Function that creates the needed high/low columns
create_high_low <- function(i, threshold, var_of_interest, grouping, 
                            dataset = stool_data){
  # i is the study
  # threshold is the vector of median values for alpha measures of interest
  # var_of_interest is the alpha metrics being used
  # grouping is the name of the case/control column
  # dataset is default to the stool_data list to allow for mapply to work
  
  # get specific data table of interest based on study
  select_data <- dataset[[i]]
  
  # create a vector with high/low versus the median value provided
  high_low <- factor(ifelse(select_data[, var_of_interest] <= threshold, 
                     invisible("low"), invisible("high")), levels = c("low", "high"))
  # Returns the vector of high/low calls
  return(high_low)
}


# Function that runs relative risk test on single variable
run_rr <- function(high_low_vector, disease_vector){
  # high_low_vector is the respective call columns from high_low for a specific alpha measure
  # disease_vector is the "is_cancer" vector is case/control info
  
  # Creates a 2x2 table of counts
  contingency <- table(high_low_vector, disease_vector)
  # runs the RR test based on the obtained 2x2 table
  test <- epi.2by2(contingency, method="cohort.count")
  # Pull only specific information from the stored list in "test"
  test_values <- cbind(test$massoc$OR.strata.score, 
                       pvalue = test$massoc$chisq.strata$p.value)
  # store both the obtained raw counts and the resulting RR with pvalue
  combined_data <- list(data_tbl = contingency, test_values = test_values)
  # Returns a list with all information needed for downstream analysis
  return(combined_data)
}


# Function to seperate out the table data from the individual analysis data
pull_data <- function(var_of_int, i, result, datalist = ind_study_data){
  # var_of_int is the alpha measures used e.g. "sobs"
  # i is the study
  # result is the type of data we want either "test_values" or "data_tbl"
  # datalist is defaulted to ind_study_data to make it easier to work with mapply
  
  # Pull the needed data and add identifiers
  tempData <- datalist[[i]][[var_of_int]][[result]] %>% as.data.frame() %>% 
    mutate(measure = var_of_int, study = i)
  # return the pulled data
  return(tempData)
  
}


# A control function to direct final table creation for ind RR analysis
make_list <- function(i, result, datalist = ind_study_data){
  # i is the study
  # result is they type of data being pulled "test_values" or "data_tbl"
  # datalist is defaulted to ind_study_data to make it easier to work with mapply
  
  # runs the function iteratively to collect the specific data
  pulled_data <- mapply(pull_data, c("sobs", "shannon", "shannoneven"), 
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
  rr_pooled_test <- rma(ai = low_Y, bi = low_N, 
                        ci = high_Y, di = high_N, data = test_data, 
                        measure = "OR", method = "REML")
  # Store a vector of the important results of interest
  results <- c(exp(c(rr = rr_pooled_test$b[[1, 1]], ci_lb = rr_pooled_test$ci.lb, 
               ci_ub=rr_pooled_test$ci.ub)), pvalue = rr_pooled_test$pval, 
               measure = alpha_d)
  # returns the vector of results
  return(results)
  
}

###########################################################################################
######### Execute the functions and munge data to "pretty" format as needed ###############
###########################################################################################


# Read in the respective data
stool_data <- mapply(get_data, c(stool_sets, both_sets), "stool", SIMPLIFY = F)

# Generate RR and data tables for every study
ind_study_data <- mapply(analyze_study, c(stool_sets, both_sets), "disease", SIMPLIFY = F)

# Pull out the RR for every study
ind_RR_data <- mapply(make_list, c(stool_sets, both_sets), "test_values", SIMPLIFY = F) %>% 
  bind_rows()

# Pull out the counts for every study and 
# merge the two different grouping columns together (is.cancer Y/N and high_low low/high)
ind_counts_data <- mapply(make_list, c(stool_sets, both_sets), "data_tbl", SIMPLIFY = F) %>% 
  bind_rows() %>% unite(group, high_low_vector, disease_vector, sep = "_") %>% 
  spread(group, Freq)

# Run pooled test
pooled_results <- t(mapply(run_pooled, c("sobs", "shannon", "shannoneven"), USE.NAMES = F)) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  mutate_at(c("rr", "ci_lb", "ci_ub", "pvalue"), as.numeric)
  

# Write out the important tables
write.csv(ind_counts_data, "data/process/tables/alpha_group_counts_summary.csv", row.names = F)
write.csv(ind_RR_data, "data/process/tables/alpha_OR_ind_results.csv", row.names = F)
write.csv(pooled_results, "data/process/tables/alpha_OR_composite.csv", row.names = F)


