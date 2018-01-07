### Run Random Forest Analysis OTU Level -- select OTUS
### Generate model and then test on remaining studies
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret", "pROC"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
# Ignore weir since none have nzv
stool_sets <- c("wang", "ahn", "zeller", "baxter", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")

# Specific Genera OTUs should belong to
rr_data <- read_csv("data/process/tables/select_genus_RR_stool_composite.csv") %>% arrange(pvalue, rr)

top5_pos_RR <- as.data.frame(rr_data %>% filter(rr > 1) %>% slice(1:5) %>% select(measure))[, "measure"]
top5_neg_RR <- as.data.frame(rr_data %>% filter(rr < 1) %>% slice(1:5) %>% select(measure))[, "measure"]
crc_genera <- c(top5_pos_RR, top5_neg_RR)

# Function to read in taxonomies and pull specific OTUs within specific Taxa
generate_select_OTUS <- function(study, specific_genera, file_path, ending){
  
  tempData <- read_tsv(paste(file_path, study, "/", study, ending, sep = "")) %>% 
    mutate(Taxonomy = gsub("\\(\\d*\\)", "", Taxonomy)) %>% 
    separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
    filter(Genus %in% specific_genera)
  
  tempList <- as.data.frame(tempData)[, "OTU"]
  
  
  return(tempList)
  
}

# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_data <- function(i, otu_list){
  # i is the study of interest
  
  # grabs subsampled data and assigns rownames from sample names to table
  shared_data <- read_tsv(paste("data/process/", i, "/", i, ".0.03.subsample.shared", 
                                sep = "")) %>% 
    select(-label, -numOtus) %>% 
    select(Group, one_of(otu_list[[i]]))
  
  # grabs the meta data and transforms polyp to control (polyp/control vs cancer) 
  study_meta <- get_file(i, "data/process/", ".metadata", rows_present = F,  
                         "stool", metadata = T) %>% filter(disease != "polyp")
  
  # Looks for Na in the meta data of interest and removes respective samples
  study_meta <- study_meta %>% filter(!is.na(disease))
  
  # conditional that checks for whether length of rows of meta data is smaller
  if(length(rownames(study_meta)) < length(rownames(shared_data))){
    # grab only the samples in the meta data file for down stream analysis
    shared_data <- shared_data %>% slice(match(study_meta$sampleID, Group))
    
    if(i == "hale"){
      
      study_meta <- study_meta %>% slice(match(shared_data$Group, sampleID))
    }
    
  } else{
    # grab only files in the data file for analysis
    study_meta <- study_meta %>% slice(match(shared_data$Group, sampleID))
  }
  # Prints out the total number of genera for that specific study
  print(paste("Total number of columns in", i, "is", 
              length(colnames(shared_data))))
  # creates a list file with both data sets
  dataList <- list(shared_data = shared_data, 
                   study_meta = study_meta, 
                   column_length = length(colnames(shared_data)))
  # returns the combined list file
  return(dataList)
  
}


# Function that grabs the meta data and replaces sampleID with disease call
assign_disease <- function(metadata_table_name, 
                           shared_data_name, fullDataList, randomize = "include"){
  # metadata_table_name is the variable with the name of the metadata file
  # shared_data_name is the variable with the name of the shared file
  # fullDataList is the original created data list
  
  # Get the respective metadata file of interest
  tempMetadata <- fullDataList[[metadata_table_name]]
  
  # create a random group label
  vars_to_sample <-  ifelse(tempMetadata$disease != "cancer", invisible(0), invisible(1))
  set.seed(12345)
  random_sample <- sample(vars_to_sample)
  
  
  # Gets transforms sample_ID column into a disease column with control v cancer calls
  tempData <- fullDataList[[shared_data_name]] %>% 
    mutate(Group = factor(ifelse(tempMetadata$disease == "polyp", 
                                 invisible("control"), 
                                 ifelse(tempMetadata$disease == "normal", 
                                        invisible("control"), 
                                        invisible(tempMetadata$disease))), 
                          levels = c("control", "cancer")), 
           random_disease = factor(ifelse(random_sample == 1, invisible("cancer"), 
                                          invisible("control")), 
                                   levels = c("control", "cancer"))) %>% 
    rename(disease = Group) %>% select(disease, random_disease, everything())
  # Returns the modified data frame that can be used for RF analysis
  return(as.data.frame(tempData))
  
}


# Function to apply and get the nzv and preProcess for the training data
get_align_info <- function(datatable){
  # datatable is the RF data table (OTU + disease + random) for study of interest
  
  # stores the disease vector (it gets removed during processing for some studies)
  disease <- datatable$disease
  random_disease <- datatable$random_disease
  # gets the respective data set i for training
  training_data <- datatable %>% select(-disease, -random_disease)
  # Check for columns that have near zero variance
  nzv <- nearZeroVar(training_data)
  # check to see if at least one value has near zero variance
  if(length(nzv) == 0 | length(nzv) == length(colnames(training_data)) | 
     length(nzv) == length(colnames(training_data)) - 1){
    # No nzv then assign training data to be itself
    training_data <- training_data
  } else{
    
    # remove columns that have near zero variance
    training_data <- training_data[, -nzv]
  }
  
  # Re add disease to the training data at the beginning of the data table
  train_data <- training_data %>% 
    mutate(disease = disease) %>% 
    select(disease, everything())
  # Re add random_disease to the random data at the beginning of the data table
  random_data <- training_data %>% 
    mutate(disease = random_disease) %>% 
    select(disease, everything())
  # create a final list with the tranformed data, the nzv columns, and the transformations
  final_info <- list(train_data = train_data, 
                     rand_data = random_data)
  # Write out the final data list
  return(final_info)
}


# Function that will run and create the needed models
make_rf_model <- function(run_marker, study, train_data){
  # run_marker is the model iteration that has currently completed
  # study is a vector for the data set to be used
  # train_data is the data table to be used for model training
  
  # Create a smaller cross-validation amount based on small n present
  if(study %in% c("weir")){
    # type of training method to be used (cv stands for cross-validation)
    method_used <- "cv"
    
    #Create Overall specifications for model tuning
    # number controls fold of cross validation
    # Repeats control the number of times to run it
    
    fitControl <- trainControl(## 10-fold CV
      method = method_used,
      number = 2,
      p = 0.8,
      classProbs = TRUE, 
      summaryFunction = twoClassSummary, 
      savePredictions = "final")
    # Uses a default 10 fold cross validation if n is not small  
  } else{
    # type of method to be used (cv stands for cross-validation)
    method_used <- "cv"
    
    #Create Overall specifications for model tuning
    # number controls fold of cross validation
    # Repeats control the number of times to run it
    
    fitControl <- trainControl(## 10-fold CV
      method = method_used,
      number = 10,
      p = 0.8, 
      classProbs = TRUE, 
      summaryFunction = twoClassSummary, 
      savePredictions = "final")
  }
  
  
  
  # Set the mtry to be based on the number of total variables in data table to be modeled
  # this formula seems to be an accepted default to use
  number_try <- round(sqrt(ncol(select(train_data, contains("Otu")))))
  
  # Set the mtry hyperparameter for the training model
  tunegrid <- expand.grid(.mtry = number_try)
  
  #Train the model
  training_model <- 
    train(disease ~ ., data = train_data, 
          method = "rf", 
          ntree = 500, 
          trControl = fitControl,
          tuneGrid = tunegrid, 
          metric = "ROC", 
          na.action = na.omit, 
          verbose = FALSE)
  
  #Print out tracking message
  print(paste("Completed ", run_marker, " RF model for ", 
              study, " using ", method_used, sep = ""))
  
  # Return the model object
  return(training_model)
}


# Function to get the min and max models to generate roc curves for
# actual and random models.
get_min_max <- function(a_models, r_models, a_summary, r_summary){
  # a_models is an object with the actual models
  # r_models is an object with the random models
  # a_summary is a table with the actual model information for each run (e.g. AUC, sens, etc.)
  # r_summary is a table with the random model information for each run (e.g. AUC, sens, etc.)
  
  # get the min and max AUC for the actual model
  a_min_row <- as.numeric((a_summary %>% filter(ROC == min(ROC)) %>% select(runs))[, "runs"])
  a_max_row <- as.numeric((a_summary %>% filter(ROC == max(ROC)) %>% select(runs))[, "runs"])
  # Get the min and max AUC for the random model
  r_min_row <- as.numeric((r_summary %>% filter(ROC == min(ROC)) %>% select(runs))[, "runs"])
  r_max_row <- as.numeric((r_summary %>% filter(ROC == max(ROC)) %>% select(runs))[, "runs"])
  # Check to see if there are more than one option for best or worse AUC
  if(length(r_max_row) > 1 | length(r_min_row) > 1 | 
     length(a_max_row) > 1 | length(a_min_row) > 1){
    # Takes the first choice if there are multiple to choose from
    a_min_row <- a_min_row[1]
    a_max_row <- a_max_row[1]
    r_min_row <- r_min_row[1]
    r_max_row <- r_max_row[1]
    
  }
  
  # Create a summary list with the best and worst model information 
  # for the actual and random models
  tempList <- list(
    actual_mod = list(
      min_model = a_models[[a_min_row]], 
      max_model = a_models[[a_max_row]]), 
    random_mod = list(
      min_model = r_models[[r_min_row]], 
      max_model = r_models[[r_max_row]]))
  # Return the variable to the global work environment
  return(tempList)
}



# Function that generates ROC curves and then compares them to random
make_summary_data <- function(i, model_info, dataList, a_summary, r_summary,  
                              train_name, random_name){
  # i is the study variable
  # model_info is a list with the best and worst models for the actual and random models
  # dataList is the list that contains the acutal data that was used to create the models
  # a_summary is a table with the actual model information for each run (e.g. AUC, sens, etc.)
  # r_summary is a table with the random model information for each run (e.g. AUC, sens, etc.)
  # train_name is the name of the training set
  # random_name is the name of the random set
  
  # Generate the best and worst roc curves for the actual model
  best_actual_roc <- roc(dataList[[train_name]]$disease ~ 
                           model_info[["actual_mod"]][["max_model"]][["pred"]][, "cancer"])
  worst_actual_roc <- roc(dataList[[train_name]]$disease ~ 
                            model_info[["actual_mod"]][["min_model"]][["pred"]][, "cancer"])
  
  # Generate the best and worst roc curves for the random model
  best_random_roc <- roc(dataList[[random_name]]$disease ~ 
                           model_info[["random_mod"]][["max_model"]][["pred"]][, "cancer"])
  worst_random_roc <- roc(dataList[[random_name]]$disease ~ 
                            model_info[["random_mod"]][["min_model"]][["pred"]][, "cancer"])
  # Generate a p-value on whether the distribution between actual and random are different
  if(i != "hale"){
    
    pvalue <- t.test(a_summary$ROC, r_summary$ROC)$p.value
  } else{
    
    pvalue <- NA
  }
  
    # Create a final list with all the needed ROC curve data and respective p-value
  finalData <- list(
    all_data = cbind(
      sens = c(best_actual_roc$sensitivities, worst_actual_roc$sensitivities, 
               best_random_roc$sensitivities, worst_random_roc$sensitivities), 
      spec = c(best_actual_roc$specificities, worst_actual_roc$specificities, 
               best_random_roc$specificities, worst_random_roc$specificities), 
      type = c(rep("actual_mod", 
                   length(c(best_actual_roc$sensitivities, worst_actual_roc$sensitivities))), 
               rep("random_mod", 
                   length(c(best_random_roc$sensitivities, worst_random_roc$sensitivities)))), 
      roc_type = c(rep("best", length(best_actual_roc$sensitivities)), 
                   rep("worst", length(worst_actual_roc$sensitivities)), 
                   rep("best", length(best_random_roc$sensitivities)), 
                   rep("worst", length(worst_random_roc$sensitivities)))) %>% 
      as.data.frame(., stringsAsFactors = F) %>% 
      mutate(sens = as.numeric(sens), spec = as.numeric(spec), 
             study = rep(i, length(spec))), 
    pvalue = pvalue)
  # Return the final list back to the global work environment
  return(finalData)
  
}



##############################################################################################
############### Run the actual programs to get the data (ALL Data) ###########################
##############################################################################################

# Set up storage variables
all_roc_data <- NULL
all_comparisons <- NULL

# Set up direction variables to set number of models to run
actual_runs <- paste("act_model_", seq(1:100), sep = "")
random_runs <- paste("rand_model_", seq(1:100), sep = "")

select_OTUs <- sapply(c(stool_sets, "flemer"), function(x) generate_select_OTUS(
  x, crc_genera, "data/process/", ".taxonomy"), simplify = F)

# Iteratively run through each study for stool
for(i in c(stool_sets, "flemer")){
  # Gets the respective data
  dataList <- get_data(i = i, select_OTUs)
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
write.csv(all_roc_data, "data/process/tables/stool_rf_select_otu_roc.csv", row.names = F)
write.csv(all_comparisons, "data/process/tables/stool_rf_select_otu_random_comparison_summary.csv", 
          row.names = F)















