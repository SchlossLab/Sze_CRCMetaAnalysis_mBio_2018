### Get Stool relative risk comparisons
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "epiR", "metafor"))


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


# Function to read in the respective transformed data tables for each study
# These should have been previously power transformed to correct skew
get_data <- function(i, sampleType){
  # i represents the data set
  # sampleType represents whether it is stool or tissue
  
  # Command that actually does the reading in
  data_list <- read.csv(paste("data/process/tables/", i, "_", sampleType, "_",  
                              "alpha_raw_values.csv", sep = ""), 
                        header = T, stringsAsFactors = F)
  
  # return to working environment the data list
  return(data_list)
}


# Analyze the data with respect to high low column table
analyze_study <- function(i, group_column, dataset = stool_data){
  
  is_cancer <- factor(ifelse(dataset[[i]][, group_column] == "cancer", 
                      invisible("Y"), invisible("N")), levels = c("Y", "N"))
  
  thresholds <- apply(select(dataset[[i]], one_of("sobs", "shannon", "shannoneven")), 2, 
                      function(x) median(x))
  
  highs_lows <- mapply(create_high_low, i, thresholds, 
                       c("sobs", "shannon", "shannoneven"), "disease", SIMPLIFY = F)
  names(highs_lows) <- c("sobs", "shannon", "shannoneven")
  
  obtained_rr <- lapply(highs_lows, 
                        function(x) run_rr(high_low_vector = x, disease_vector = is_cancer)) 
  return(obtained_rr)
}


# create the needed high/low columns
create_high_low <- function(i, threshold, var_of_interest, grouping, 
                            dataset = stool_data){
  
  select_data <- dataset[[i]]
  
  high_low <- factor(ifelse(select_data[, var_of_interest] <= threshold, 
                     invisible("low"), invisible("high")), levels = c("low", "high"))
  
  return(high_low)
}


# Run relative risk test on single variable
run_rr <- function(high_low_vector, disease_vector){
  
  contingency <- table(high_low_vector, disease_vector)
  
  test <- epi.2by2(contingency, method="cohort.count")
  
  test_values <- cbind(test$massoc$RR.strata.score, 
                       pvalue = test$massoc$chisq.strata$p.value)
  
  combined_data <- list(data_tbl = contingency, test_values = test_values)
  
  return(combined_data)
}


# Seperate out the table data from the individual analysis data
pull_data <- function(var_of_int, i, result, datalist = ind_study_data){
  
  tempData <- datalist[[i]][[var_of_int]][[result]] %>% as.data.frame() %>% 
    mutate(measure = var_of_int, study = i)
  
  return(tempData)
  
}


# Control function to direct table creation
make_list <- function(i, result, datalist = ind_study_data){
  
  pulled_data <- mapply(pull_data, c("sobs", "shannon", "shannoneven"), 
                 i, result, SIMPLIFY = F) %>% bind_rows()
  
  return(pulled_data)
  
  
}


# Function to run test for selected alpha measure
run_pooled <- function(alpha_d, dataset = ind_counts_data){
  
  test_data <- dataset %>% filter(measure == alpha_d)
  
  rr_pooled_test <- rma(ai = low_Y, bi = low_N, 
                        ci = high_Y, di = high_N, data = test_data, 
                        measure = "RR", method = "REML")
  
  results <- c(exp(c(rr = rr_pooled_test$b[[1, 1]], ci_lb = rr_pooled_test$ci.lb, 
               ci_ub=rr_pooled_test$ci.ub)), pvalue = rr_pooled_test$pval, 
               measure = alpha_d)
  
  return(results)
  
}



# Read in the respective data
stool_data <- mapply(get_data, c(stool_sets, both_sets), "stool", SIMPLIFY = F)

# Generate RR and data tables for every study
ind_study_data <- mapply(analyze_study, c(stool_sets, both_sets), "disease", SIMPLIFY = F)

# Pull out the RR for every study
ind_RR_data <- mapply(make_list, c(stool_sets, both_sets), "test_values", SIMPLIFY = F) %>% 
  bind_rows()

# Pull out the counts for every study
ind_counts_data <- mapply(make_list, c(stool_sets, both_sets), "data_tbl", SIMPLIFY = F) %>% 
  bind_rows() %>% unite(group, high_low_vector, disease_vector, sep = "_") %>% 
  spread(group, Freq)

# Run pooled test
pooled_results <- t(mapply(run_pooled, c("sobs", "shannon", "shannoneven"), USE.NAMES = F)) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  mutate_at(c("rr", "ci_lb", "ci_ub", "pvalue"), as.numeric)
  

# Write out the important tables
write.csv(ind_counts_data, "data/process/tables/alpha_group_counts_summary.csv", row.names = F)
write.csv(ind_RR_data, "data/process/tables/alpha_RR_ind_results.csv", row.names = F)
write.csv(pooled_results, "data/process/tables/alpha_RR_composite.csv", row.names = F)










###### TO DO LIST ######

## Need to seperate the respective tables from the results 
## Need to format the table results to be used in an aggregate analysis
## Run aggregate analysis on every measure


