### Get Tissue alpha comparisons
### One for matched and one for unmatched
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "lme4"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

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


# Run tests for alpha metrics using t-tests
get_ttest_comparisons <- function(alpha_metric, split_on, data_set, paired){
  # alpha_metric represents the alpha variables of interest
  # split_on represents how to split the data for the t-test
  # data_set is the new combined data table (created separately)
  # coded like this to make it easier to use with mapply
  
  data_set <- get(data_set)
  
  # Performs the actual test
  temp_test <- t.test(
    filter(data_set, disease == split_on)[, alpha_metric], 
    filter(data_set, disease != split_on)[, alpha_metric],
                      alternative = "two.sided", 
                      paired = ifelse(paired == "Yes", TRUE, FALSE))
  
  # This creates a easy to understand table of the results
  test_info <- data_frame(c(temp_test$estimate, temp_test$p.value)) %>% 
    setNames("results") %>% 
    mutate(col_names = if(paired == "Yes"){
      c("mean_difference", "pvalue")}else{
        c(split_on, paste("non_", split_on, sep = ""), "pvalue")
      }) %>% 
    spread(col_names, results) %>% mutate(measures = alpha_metric)
  
  # Return the easy to understand table
  return(test_info)
  
}


#Run an ANOVA with Tukey post hoc test to look for differences within each group
get_anova_comparisons <- function(alpha_metric, set_groups = "disease", 
                                  data_set = tissue_unmatched){
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
                             individual = NULL, data_set){
  # alpha_metric represents the alpha variables of interest
  # study_group represents the column indicating which study the data belongs to
  # disease_group represents the column where the group variables are listed
  # data_set is the combined data set (used in all testing functions)
  
  data_set <- get(data_set)
  
  # Create the null model without group of interest
  null_model <- lmer(
    if(is.null(individual)){
      
      as.formula(paste(alpha_metric, " ~ ", "(1|", study_group, ")", sep = ""))
    } else{
      
      as.formula(paste(alpha_metric, " ~ ", "(1|", study_group, ") + ", 
                       "(1|", individual, ")", sep = ""))
    }, data = data_set, REML=FALSE)
  
  # create the disease model with the group of interest
  disease_model <- lmer(
    if(is.null(individual)){
      
      as.formula(paste(alpha_metric, " ~ ", "(1|", disease_group, ") + ", 
                       "(1|", study_group, ")", sep = ""))
    } else{
      
      as.formula(paste(alpha_metric, " ~ ", "(1|", disease_group, ") + ", 
                       "(1|", study_group, ") + ", "(1|", individual, ")", 
                       sep = ""))
    }, data = data_set, REML=FALSE)
  
  # Test whether the disease model explains significantly more variation than the null model
  results <- anova(null_model, disease_model)
  
  # Pull out metrics of interests
  imp_values <- c(chi_sq = results[, "Chisq"][2], pvalue = results[, "Pr(>Chisq)"][2])
  
  # Return to environment the metrics of interst
  return(imp_values)
  
}


# Function to write out all the results
make_the_tables <- function(table_name){
  # table_name is the name of the table to get to write out
  
  data_table <- get(table_name)
  
  # saves the data table outside of R
  write.csv(data_table, paste("data/process/tables/alpha_", table_name, ".csv"), 
            row.names = F)
}



# Filter out data that are for polyp and control matched samples
samples_to_remove <- tissue_matched %>% 
  filter(disease == "polyp" | disease == "Polyp") %>% 
  select(id) %>% 
  bind_rows(tissue_matched %>% filter(disease == "control") %>% select(id, matchings) %>% 
    group_by(id) %>% filter(n() > 1) %>% unique() %>% select(id))

good_tissue_matched <- tissue_matched %>% 
  filter(!as.character(id) %in% as.character(samples_to_remove[, "id"]))


# Run ttest
unmatched_ttest_tissue <- bind_rows(mapply(
  get_ttest_comparisons, c("sobs", "shannon", "shannoneven"), 
  "cancer", "tissue_unmatched", "No", USE.NAMES = T, SIMPLIFY = F)) %>% 
  mutate(BH = p.adjust(pvalue, method = "BH")) %>% 
  select(cancer, non_cancer, pvalue, BH, measures)

matched_ttest_tissue <- bind_rows(mapply(
  get_ttest_comparisons, c("sobs", "shannon", "shannoneven"), 
  "cancer", "good_tissue_matched", "Yes", USE.NAMES = T, SIMPLIFY = F)) %>% 
  mutate(BH = p.adjust(pvalue, method = "BH")) %>% 
  select(mean_difference, pvalue, BH, measures)


# Run ANOVA with Tukey
unmatched_tukey_results <- mapply(get_anova_comparisons, c("sobs", "shannon", "shannoneven"), 
                        USE.NAMES = T, SIMPLIFY = F) %>% bind_rows()


# RUn linear mixed-effect models
unmatched_mixeffect_results <- mapply(get_mixed_effect, c("sobs", "shannon", "shannoneven"), 
                            data_set = "tissue_unmatched", USE.NAMES = T, SIMPLIFY = F) %>% 
  bind_rows() %>% gather(key = alpha_metric, value = value) %>% 
  mutate(measure = rep(c("chi_sq", "pvalue"), length(value)/2)) %>% 
  spread(measure, value) %>% 
  mutate(BH = p.adjust(pvalue, method = "BH"))


matched_mixeffect_results <- mapply(get_mixed_effect, c("sobs", "shannon", "shannoneven"), 
                                      individual = "id", data_set = "good_tissue_matched", 
                                    USE.NAMES = T, SIMPLIFY = F) %>% 
  bind_rows() %>% gather(key = alpha_metric, value = value) %>% 
  mutate(measure = rep(c("chi_sq", "pvalue"), length(value)/2)) %>% 
  spread(measure, value) %>% 
  mutate(BH = p.adjust(pvalue, method = "BH"))


# Write out the results from the 3 tests
mapply(make_the_tables, 
       c("unmatched_ttest_tissue", "matched_ttest_tissue", "unmatched_tukey_results", 
         "unmatched_mixeffect_results", "matched_mixeffect_results"))





