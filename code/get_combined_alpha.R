### Generate Alpha Diversity Comparisons
### Z transform power transformed data
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "lme4"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("wang", "brim", "weir", "ahn", "zeller", "baxter")

# Both Tissue and Stool
  # flemer sampletype = biopsy or stool
  # chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")


# Function to read in the respective transformed data tables for each study
get_transformed_data <- function(i, sampleType){
  
  data_list <- read.csv(paste("data/process/tables/", i, "_", sampleType, "_",  
                                   "transformed_data_alpha_raw_values.csv", sep = ""), 
                             header = T, stringsAsFactors = F)

  return(data_list)
}


# Function to apply a z-score normalization and and a new column with these values
# for sobs, shannon, and shannoneven
zscore_transform <- function(dataList, 
                             variables = c("sobs", "shannon", "shannoneven")){
  
  z_trans_data <- dataList %>% 
    mutate_at(vars(variables), funs(as.numeric(scale(.))))
  
  return(z_trans_data)
}


# Run tests for alpha metrics using t-tests
get_ttest_comparisons <- function(alpha_metric, split_on, data_set = combined_data){
  
  temp_test <- t.test(filter(data_set, disease == split_on)[, alpha_metric], 
         filter(data_set, disease != split_on)[, alpha_metric], 
         alternative = "two.sided")

  test_info <- data_frame(c(temp_test$estimate, temp_test$p.value)) %>% 
    setNames("results") %>% 
    mutate(col_names = c(split_on, paste("non_", split_on, sep = ""), "pvalue")) %>% 
    spread(col_names, results) %>% mutate(measures = alpha_metric)

  return(test_info)
  
}

#Run an ANOVA with Tukey post hoc test to look for differences within each group
get_anova_comparisons <- function(alpha_metric, set_groups = "disease", 
                                  data_set = combined_data){
  
    test_data <- TukeyHSD(aov(lm(
    as.formula(paste(alpha_metric, "~", set_groups)), data = data_set)))[[set_groups]] %>% 
    as.data.frame() %>% 
    mutate(comparison = rownames(.), measure = alpha_metric)
  
  return(test_data)
  
}


# Run a linear mixed-effect model accounting for study
# for each alpha metric used
get_mixed_effect <- function(alpha_metric, study_group = "study", disease_group = "disease", 
                             data_set = combined_data){
  
  null_model <- lmer(
    as.formula(paste(alpha_metric, " ~ ", "(1|", study_group, ")", sep = "")), 
    data = data_set, REML=FALSE)
  
  disease_model <- lmer(
    as.formula(paste(alpha_metric, " ~ ", "(1|", disease_group, ") + ", 
                     "(1|", study_group, ")", sep = "")), 
    data = data_set, REML=FALSE)
  
  results <- anova(null_model, disease_model)
  
  imp_values <- c(chi_sq = results[, "Chisq"][2], pvalue = results[, "Pr(>Chisq)"][2])
  
  return(imp_values)
  
}


# Read in the respective data
stool_transformed_data <- mapply(get_transformed_data, 
                                 c(stool_sets, both_sets), "stool", SIMPLIFY = F)

# Apply z-score normalization to every data set
stool_ztrans_pwrtrans_data <- lapply(stool_transformed_data, zscore_transform)

# create a combined data table
combined_data <- bind_rows(lapply(stool_ztrans_pwrtrans_data, 
                                  function(x) as.data.frame(x) %>% 
                           mutate(group = as.character(group))))

### Need to loop this for each respective variable ###

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
make_the_tables <- function(i, table_name){
  
  write.csv(i, paste("data/process/tables/alpha_", table_name, ".csv"), row.names = F)
}

mapply(make_the_tables, 
       list(ttest_results, tukey_results, mixeffect_results), 
       c("ttest_results", "tukey_results", "mixeffect_results"))













