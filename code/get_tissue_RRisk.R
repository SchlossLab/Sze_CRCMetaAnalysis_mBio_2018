### Get Tissue relative risk comparisons
### One for unmatched and one for matched if possible
### Marc Sze

# For this analysis combined both unmatched and matched together to boost n and power
# Can label which studies had majority matched to see if that made a difference at all.

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


# Read in specific data tables to be used

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0))

tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  mutate(disease = gsub("adenoma", "polyp", disease))


combined_tissue <- tissue_unmatched %>% 
  select(one_of("group", "study", "disease", "r_sobs", "r_shannon", "r_shannoneven")) %>% 
           bind_rows(
             select(tissue_matched, 
                    one_of("group", "study", "disease", "r_sobs", "r_shannon", "r_shannoneven")))

# Remove polyp only group
no_p_tissue_unmatched <- combined_tissue %>% filter(study != "lu")

# Function to run the analysis
analyze_study <- function(i, group_column, dataset = no_p_tissue_unmatched){
  
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
  # Obtains the individual relative risk and CI for each study
  obtained_rr <- lapply(highs_lows, 
                        function(x) run_rr(high_low_vector = x, disease_vector = is_cancer))
  
  return(obtained_rr)
  
}


# Function that creates the needed high/low columns
create_high_low <- function(i, var_of_interest, threshold, grouping, 
                            dataset = no_p_tissue_unmatched){
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


# Function that runs relative risk test on single variable
run_rr <- function(high_low_vector, disease_vector){
  # high_low_vector is the respective call columns from high_low for a specific alpha measure
  # disease_vector is the "is_cancer" vector is case/control info
  
  # Creates a 2x2 table of counts
  contingency <- table(high_low_vector, disease_vector)
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


# Generate RR and data tables for every study
ind_study_data <- mapply(analyze_study, c(tissue_sets, both_sets), "disease", SIMPLIFY = F)

# Pull out the RR for every study
ind_RR_data <- mapply(make_list, c(tissue_sets, both_sets), "test_values", SIMPLIFY = F) %>% 
  bind_rows()


# Pull out the counts for every study and 
# merge the two different grouping columns together (is.cancer Y/N and high_low low/high)
ind_counts_data <- mapply(make_list, c(tissue_sets, both_sets), "data_tbl", SIMPLIFY = F) %>% 
  bind_rows() %>% unite(group, high_low_vector, disease_vector, sep = "_") %>% 
  spread(group, Freq)


# Run pooled test
pooled_results <- t(mapply(run_pooled, c("sobs", "shannon", "shannoneven"), USE.NAMES = F)) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  mutate_at(c("rr", "ci_lb", "ci_ub", "pvalue"), as.numeric)


# Write out the important tables
write.csv(ind_counts_data, "data/process/tables/alpha_group_counts_tissue_summary.csv", row.names = F)
write.csv(ind_RR_data, "data/process/tables/alpha_RR_ind_tissue_results.csv", row.names = F)
write.csv(pooled_results, "data/process/tables/alpha_RR_tissue_composite.csv", row.names = F)



