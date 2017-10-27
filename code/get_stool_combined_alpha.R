### Generate Alpha Diversity Comparisons
### Z transform power transformed data
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "car", "lme4"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("wang", "brim", "weir", "ahn", "zeller", "baxter", "hale")

# Both Tissue and Stool
  # flemer sampletype = biopsy or stool
  # chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")


# Function to read in the respective transformed data tables for each study
# These should have been previously power transformed to correct skew
get_transformed_data <- function(i, sampleType){
  # i represents the data set
  # sampleType represents whether it is stool or tissue
  
  # Command that actually does the reading in
  data_list <- read.csv(paste("data/process/tables/", i, "_", sampleType, "_",  
                                   "transformed_data_alpha_raw_values.csv", sep = ""), 
                             header = T, stringsAsFactors = F)
  
  # return to working environment the data list
  return(data_list)
}


# Function to apply a z-score normalization and and a new column with these values
# for sobs, shannon, and shannoneven
zscore_transform <- function(dataList, 
                             variables = c("sobs", "shannon", "shannoneven")){
  # dataList represents the list data that was read in by get_transformed_data
  # variables reprents data of interest for downstream comparisons
  
  # command that z-score transforms interested data 
  z_trans_data <- dataList %>% 
    mutate_at(vars(variables), funs(as.numeric(scale(.))))
  
  # return to working environment the new z-score transformed data
  return(z_trans_data)
}


# Run tests for alpha metrics using t-tests
get_ttest_comparisons <- function(alpha_metric, split_on, data_set = combined_data){
  # alpha_metric represents the alpha variables of interest
  # split_on represents how to split the data for the t-test
  # data_set is the new combined data table (created separately)
    # coded like this to make it easier to use with mapply
  
  # Performs the actual test
  temp_test <- t.test(filter(data_set, disease == split_on)[, alpha_metric], 
         filter(data_set, disease != split_on)[, alpha_metric], 
         alternative = "two.sided")
  
  # This creates a easy to understand table of the results
  test_info <- data_frame(c(temp_test$estimate, temp_test$p.value)) %>% 
    setNames("results") %>% 
    mutate(col_names = c(split_on, paste("non_", split_on, sep = ""), "pvalue")) %>% 
    spread(col_names, results) %>% mutate(measures = alpha_metric)
  
  # Return the easy to understand table
  return(test_info)
  
}

#Run an ANOVA with Tukey post hoc test to look for differences within each group
get_anova_comparisons <- function(alpha_metric, set_groups = "disease", 
                                  data_set = combined_data){
  # alpha_metric represents the alpha variables of interest
  # set_groups specifies what column to use for the groups variable
  # data_set is the combined data set (used in all testing functions)
    
  # this runs an ANOVA with Tukey post hoc test and makes a readable table 
  test_data <- TukeyHSD(aov(lm(
    as.formula(paste(alpha_metric, "~", set_groups)), data = data_set)))[[set_groups]] %>% 
    as.data.frame() %>% 
    mutate(comparison = rownames(.), measure = alpha_metric)
  
  # Returns the results of the test for each comparison
  return(test_data)
  
}


# Run a linear mixed-effect model accounting for study
# for each alpha metric used
get_mixed_effect <- function(alpha_metric, study_group = "study", disease_group = "disease", 
                             variable_region = "v_region", data_set = combined_data){
  # alpha_metric represents the alpha variables of interest
  # study_group represents the column indicating which study the data belongs to
  # disease_group represents the column where the group variables are listed
  # data_set is the combined data set (used in all testing functions)
  
  # Create the null model without group of interest
  null_model <- lmer(
    as.formula(paste(alpha_metric, " ~ ", "(1|", study_group, ") + ", "(1|", v_region, ")", 
                     sep = "")), data = data_set, REML=FALSE)
  
  # create the disease model with the group of interest
  disease_model <- lmer(
    as.formula(paste(alpha_metric, " ~ ", "(1|", disease_group, ") + ", 
                     "(1|", study_group, ") + ", "(1|", v_region, ")", 
                     sep = "")), data = data_set, REML=FALSE)
  
  # Test whether the disease model explains significantly more variation than the null model
  results <- anova(null_model, disease_model)
  
  # Pull out metrics of interests
  imp_values <- c(chi_sq = results[, "Chisq"][2], pvalue = results[, "Pr(>Chisq)"][2])
  
  # Return to environment the metrics of interst
  return(imp_values)
  
}


# Function to write out all the results
make_the_tables <- function(i, table_name){
  # i is the respective result data table
  # table_name is what this table should be called when saved
  
  # saves the data table outside of R
  write.csv(i, paste("data/process/tables/alpha_", table_name, ".csv", sep = ""), row.names = F)
}


#### Stool ####

# Read in the respective data
stool_transformed_data <- mapply(get_transformed_data, 
                                 c(stool_sets, both_sets), "stool", SIMPLIFY = F)

# Apply z-score normalization to every data set
stool_ztrans_pwrtrans_data <- lapply(stool_transformed_data, zscore_transform)

# create a combined data table
combined_data <- bind_rows(lapply(stool_ztrans_pwrtrans_data, 
                                  function(x) as.data.frame(x) %>% 
                           mutate(group = as.character(group)))) %>% 
  mutate(v_region = ifelse(study %in% c("baxter", "zeller", "weir"), invisible("V4"), 
                           ifelse(study %in% c("ahn", "flemer"), invisible("V3-4"), 
                                  ifelse(study %in% c("brim", "chen"), invisible("V1-V3"), 
                                         ifelse(study == "wang", invisible("V3"), invisible("V3-5"))))))


# Test the stool data for difference between cancer and non-cancer
ttest_results <- bind_rows(mapply(get_ttest_comparisons, 
                                  c("sobs", "shannon", "shannoneven"), "cancer", 
                                  USE.NAMES = T, SIMPLIFY = F)) %>% 
  mutate(BH = p.adjust(pvalue, method = "BH")) %>% 
  select(cancer, non_cancer, pvalue, BH, measures)

# Test if there is a difference based on adenoma being seperated from controls
tukey_results <- mapply(get_anova_comparisons, c("sobs", "shannon", "shannoneven"), 
               USE.NAMES = T, SIMPLIFY = F) %>% bind_rows()


# Account for data sets in analysis using a linear mixed-effect model
mixeffect_results <- mapply(get_mixed_effect, c("sobs", "shannon", "shannoneven"), 
               USE.NAMES = T, SIMPLIFY = F) %>% bind_rows() %>% 
  gather(key = alpha_metric, value = value) %>% 
  mutate(measure = rep(c("chi_sq", "pvalue"), length(value)/2)) %>% 
  spread(measure, value) %>% 
  mutate(BH = p.adjust(pvalue, method = "BH"))

# Write out the results from the 3 tests
mapply(make_the_tables, 
       list(ttest_results, tukey_results, mixeffect_results), 
       c("ttest_results", "tukey_results", "mixeffect_results"))













