### Run Random Forest Analysis OTU Level -- stool adenoma
### Generate model and then test on remaining studies
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "caret", "pROC"))


# Stool Only polyp sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("brim", "zeller", "baxter", "hale")


##############################################################################################
############### List of function to allow for the analysis to work ###########################
##############################################################################################














##############################################################################################
############### Run the actual programs to get the data (ALL Data) ###########################
##############################################################################################

# Set up storage variables
all_roc_data <- NULL
all_comparisons <- NULL

# Set up direction variables to set number of models to run
actual_runs <- paste("act_model_", seq(1:100), sep = "")
random_runs <- paste("rand_model_", seq(1:100), sep = "")

# Iteratively run through each study for stool
for(i in c(stool_sets, "flemer")){
  # Gets the respective data
  dataList <- get_data(i = i)
  # merges the needed metadata with the variables to test and creates a random label as well
  disease_dataset <- assign_disease("study_meta", "shared_data", dataList)
  # makes sure all the genera are the same for every data set to be tested
  rf_data <- get_align_info(disease_dataset)
  # generates a distribution of models based on real data
  actual_model <- sapply(actual_runs, 
                         function(x) make_rf_model(x, i, rf_data[["train_data"]]), simplify = F) 
  # generates a distribution of models based on random labeleed data
  random_model <- sapply(random_runs, 
                         function(x) make_rf_model(x, i, rf_data[["rand_data"]]), simplify = F)
  # create a data table with summary stats of all n actual models run
  actual_summary <- sapply(actual_model, 
                           function(x) x$results, simplify = F) %>% bind_rows() %>% 
    mutate(runs = rownames(.))
  # create a data table with summary stats of all n random models run
  random_summary <- sapply(random_model, 
                           function(x) x$results, simplify = F) %>% bind_rows() %>% 
    mutate(runs = rownames(.))
  # Generates the best and worst models for actual and random models
  model_info <- get_min_max(actual_model, random_model, 
                            actual_summary, random_summary)
  # Gets the pvalue between actual and random as well as needed graphing data
  test <- make_summary_data(i = i, model_info = model_info, rf_data, 
                            actual_summary, random_summary, "train_data", "rand_data")
  # transforms the collected roc data into a data frame
  all_roc_data <- all_roc_data %>% bind_rows(test[["all_data"]])
  # Creates a summary table for the comparisons made
  all_comparisons <- rbind(all_comparisons, 
                           as.data.frame.list(
                             c(actual_summary %>% summarise(act_mean_auc = mean(ROC, na.rm = T), 
                                                            act_sd_auc = sd(ROC, na.rm = T)), 
                               random_summary %>% summarise(rand_mean_auc = mean(ROC, na.rm = T), 
                                                            rand_sd_auc = sd(ROC, na.rm = T)), 
                               pvalue = test[["pvalue"]], study = i)))
  # Tracking print out that allows tracking of which study has completed
  print(paste("Completed study:", i, "RF testing"))
  
}


# Write out the relevant data frames
write.csv(all_roc_data, "data/process/tables/stool_rf_otu_roc.csv", row.names = F)
write.csv(all_comparisons, "data/process/tables/stool_rf_otu_random_comparison_summary.csv", 
          row.names = F)


