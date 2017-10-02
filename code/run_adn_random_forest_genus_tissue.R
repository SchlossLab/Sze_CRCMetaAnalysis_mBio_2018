### Run Random Forest Analysis -- adenoma tissue
### Generate model and then test on remaining studies
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "caret", "pROC"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("lu")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer")

# Read in specific data tables to be used

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0)) %>% 
  filter(study %in% c(tissue_sets, both_sets)) %>% 
  rename(sample_id = group) %>% 
  filter(id != 20) #remove the one matched control sample

tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  filter(study %in% c(tissue_sets, both_sets), disease != "cancer") %>% 
  rename(sample_id = group)


##############################################################################################
############### List of functions to run the analysis ########################################
##############################################################################################

# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_data <- function(i, metadata){
  # i is the study of interest
  
  # gets original sample names
  sample_names <- (get_file(i, "data/process/", "_genera_shared.csv") %>% 
                     mutate(sample_names = rownames(.)) %>% select(sample_names))[, "sample_names"]
  # grabs subsampled data and assigns rownames from sample names to table
  sub_genera_data <- get_file(i, "data/process/", "_subsample_genera.csv", rows_present = F, 
                              vec_of_rownames = sample_names) %>% 
    as.data.frame() %>% mutate(sample_ID = rownames(.)) %>% 
    select(sample_ID, everything()) %>% mutate(sample_ID = as.character(sample_ID))
  # grabs the meta data and transforms polyp to control (polyp/control vs cancer) 
  study_meta <- metadata %>% filter(study == i)
  
  # Looks for Na in the meta data of interest and removes respective samples
  study_meta <- study_meta %>% filter(!is.na(disease)) %>% 
    select(sample_id, disease)
  
  sub_genera_data <- study_meta %>% 
    inner_join(sub_genera_data, by = c("sample_id" = "sample_ID")) %>% 
    select(-sample_id)
  
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta, 
                   column_length = length(colnames(sub_genera_data)))
  # returns the combined list file
  return(dataList)
  
}


# Function that aligns the genera from all data sets
align_genera <- function(studies, length_column_name, 
                         genera_data_name, dataList){
  # studies is a variable for the name of study to be used
  # length_column_name is the name of the vector in dataList that stores the total genera number
  # genera_data_name is the name of the data table in dataList that contains the genera information
  # dataList is a list that has for every study the genus info, metadata, and total genera present
  
  if(length(studies) > 1){
    
    # Pulls out the total number of genera identified in each study and orders them highest to lowest
    genera_num_list <- sort.int(sapply(studies, 
                                       function(x) dataList[[x]][[length_column_name]]))
    # Counting variable (helps to direct flow)
    x = 1
    # Checks to see if the lowest value is the same as the highest, if not keep iterating through
    while(genera_num_list[1]!= genera_num_list[length(genera_num_list)]){
      # stores the  lowest genera study name
      lowest_genera_study <- names(genera_num_list[1])
      # Check to see if this is the first time through the iteration
      if(x == 1){
        # Remove any groups that have unclassified as part of their ID
        genera_names <- dataList[[lowest_genera_study]][[genera_data_name]] %>% 
          select(-contains("_unclassified")) %>% colnames(.)
      } else{
        # If not the first through get the names of genrea in the current lowest study 
        genera_names <- colnames(temp_aligned_genera[[lowest_genera_study]])
      }
      
      # iterate through each study matching only those in the lowest genera data set
      temp_aligned_genera <- suppressWarnings(sapply(studies, 
                                                     function(x) 
                                                       select(dataList[[x]][[genera_data_name]], 
                                                              one_of(genera_names)), simplify = F))
      # get the updated total genera and then sort lowest to highest
      genera_num_list <- sort.int(sapply(studies, 
                                         function(x) 
                                           length(colnames(temp_aligned_genera[[x]]))))
      # Print out an update as to how the matching is going
      print(paste("Min and Max total genera is:", 
                  min(genera_num_list), ",", max(genera_num_list)))
      # move the tracker / flow director up one
      x = x + 1
      
    }
  } else{
    
    temp_aligned_genera <- dataList[[studies]][[genera_data_name]]
    print("Only one study, simply extracting necessary data table.")
  }
  
  
  # When the while loop exits return the aligned genera data list
  return(temp_aligned_genera)
}


# Function that grabs the meta data and replaces sampleID with disease call
assign_disease <- function(studies, matched_genera){
  # studies is the variable with the names of the studies 
  # matched_genera is the data list with only genera in every study
  
  # Get the respective metadata file of interest
  tempData <- matched_genera[[studies]]
  # Gets transforms sample_ID column into a disease column with control v cancer calls
  corr_tempData <- tempData %>% 
    mutate(disease = factor(ifelse(disease == "normal", 
                                   invisible("control"), invisible(disease)), 
                            levels = c("control", "polyp")))
  # Returns the modified data frame that can be used for RF analysis
  return(corr_tempData)
  
}


# Function to apply and get the nzv and preProcess for the training data
### Errors with the pre-process so default to no preProcess for adn stool
get_align_info <- function(i, dataList){
  # i stands for the study
  # dataList is the rf_dataset set up (disease + genus info) for every study
  
  # gets the respective data set i for training
  training_data <- dataList[[i]]
  # stores the disease vector (it gets removed during processing for some studies)
  disease <- training_data$disease
  # Check for columns that have near zero variance
  nzv <- nearZeroVar(training_data)
  # remove columns that have near zero variance
  training_data <- training_data[, -nzv]
  # Check to see if disease column was removed during the processing
  if("disease" %in% colnames(training_data)){
    # keep training data the same
    training_data <- training_data
    
  } else{
    # Re add disease to the training data at the beginning of the data table
    training_data <- training_data %>% mutate(disease = disease) %>% 
      select(disease, everything())
  }
  # create a final list with the tranformed data, the nzv columns, and the transformations
  final_info <- list(train_data = training_data, 
                     near_zero_variance = nzv)
  # Write out the final data list
  return(final_info)
}


# Function to generate the preprocessing files on test data
apply_preprocess <- function(i, trainingList_info, dataList){
  # i is the study of interest
  # trainingList_info if the RF prepared training data with parameters (list)
  # dataList is all data before transformation (typically rf_datasets)
  
  # remove the trianing data from the full data list
  dataList[[i]] <- NULL
  # Pull the columns with near zero variance from the training list
  nzv <- as.numeric(trainingList_info[["near_zero_variance"]])
  # check to see if the first column is part of the nzv
  if(1 %in% nzv){
    # if 1 is part of it (this is the disease column) remove it from nzv
    nzv <- nzv[nzv != 1]
    # otherwise go here
  } else {
    # keep nzv as is 
    nzv <- nzv
  }
  # Remove near zero variance from all data sets
  dataList <- lapply(dataList, function(x) as.data.frame(x[, -nzv]))
  # return the modified list with all studies
  return(dataList)
  
}


# Function that will run and create the needed model
make_rf_model <- function(train_data){
  # train_data is the data table to be used for model training
  
  #Create Overall specifications for model tuning
  # number controls fold of cross validation
  # Repeats control the number of times to run it
  
  fitControl <- trainControl(## 10-fold CV
    method = "cv",
    number = 10,
    p = 0.8, 
    classProbs = TRUE, 
    summaryFunction = twoClassSummary, 
    savePredictions = "final")
  
  # Set the mtry to be based on the number of total variables in data table to be modeled
  # this formula seems to be an accepted default to use
  number_try <- round(sqrt(ncol(train_data)))
  
  # Set the mtry hyperparameter for the training model
  tunegrid <- expand.grid(.mtry = number_try)
  
  #Train the model
  set.seed(12345)
  training_model <- 
    train(disease ~ ., data = train_data, 
          method = "rf", 
          ntree = 500, 
          trControl = fitControl,
          tuneGrid = tunegrid, 
          metric = "ROC", 
          na.action = na.omit, 
          verbose = FALSE)
  
  # Return the model object
  return(training_model)
}


# Function that will test all existing test sets (i.e. other studies)
get_test_data <- function(i, train_study, 
                          training_model, training_data, testdataList){
  # i is the study of interest
  # train_study is the variable for the name of the study used as the training set 
  # training_model is the RF object of the training study
  # training_data is the raw data used to make the training model
  # testdataList is all the transformed data from all other studies
  
  # Get the outcomes of the voting for the training model
  train_prediction <- training_model$finalModel$votes %>% as.data.frame()
  # get the predictions for each of the data sets based on the training model
  test_predictions <- sapply(i, function(x) 
    predict(training_model, testdataList[[x]], type = 'prob'), simplify = F)
  # Generate roc curve infor (sens and spec) to be able to graph roc curves in the future
  overall_rocs <- sapply(i, function(x) 
    roc(testdataList[[x]]$disease ~ test_predictions[[x]][, "polyp"]), simplify = F)
  # add the training data roc information to this list
  overall_rocs[[train_study]] <- roc(training_data$disease ~ train_prediction[, "polyp"])
  
  # Write out all the roc information from every data set
  return(overall_rocs)
  
}

# Function that creates a finalized data table for respective model
make_data_table <- function(final_rocs){
  # final_rocs should be a list of roc objects for every study (including training)
  
  # creates a data table with relevant summary information
  tempData <- sapply(names(final_rocs), function(x) 
    as.data.frame(cbind(sens = as.numeric(final_rocs[[x]]$sensitivities), 
                        spec = final_rocs[[x]]$specificities, 
                        auc = rep(final_rocs[[x]]$auc[1], length(final_rocs[[x]]$sensitivities))), 
                  stringsAsFactors = F) %>% 
      mutate(study = rep(x, length(final_rocs[[x]]$sensitivities))), simplify = F) %>% bind_rows()
  # write out the summary information table
  return(tempData)
}



# Function to execute the major commands for RF gathering
run_rf_tests <- function(study, rf_dataList, specific_vars = F){
  # study is the study of interest
  # rf_dataList is the genera_aligned disease column added data files 
  
  if(specific_vars == F & study != "brim"){
    
    # Generate the nzv and transformations that need to be applied based on 
    # study used for training data 
    first_study <- get_align_info(study, rf_dataList)
    # remove nzv columns, transform, and normalize data sets based on training set 
    test_dataList <- apply_preprocess(study, first_study, rf_dataList)
    # Generate the RF model from the relevant training data set
    train_model_data <- make_rf_model(first_study$train_data)
    # Test the model on each of the data sets not used in training
    test_dataLists <- get_test_data(names(test_dataList), study, train_model_data, 
                                    first_study[["train_data"]], test_dataList)
    
  } else{
    
    first_study <- rf_dataList[[study]]
    test_dataList = rf_dataList
    test_dataList[[study]] <- NULL
    # Generate the RF model from the relevant training data set
    train_model_data <- make_rf_model(first_study)
    # Test the model on each of the data sets not used in training
    test_dataLists <- get_test_data(names(test_dataList), study, train_model_data, 
                                    first_study, test_dataList)
    
  }
  
  # output the final results from all tests
  return(test_dataLists)
}


# Function to generate pvalues for the results ROC cuves (default method delong)
make_model_comparisons <- function(i, rocList, comp_method = "bootstrap"){
  # i is the study of interest
  # rocList is a list of roc objects (obtained from get_test_data)
  # comp_method is a character call of what method to use for comparisons
  # the default is set to bootstrap because it can test different direction curves
  
  # create a temp list variable 
  tempList <- rocList[[i]]
  # remove the study of interest (training roc) from temp list
  tempList[[i]] <- NULL
  # create the training roc variable
  train_model_roc <- rocList[[i]][[i]]
  # iterate through each study comparing the training roc the test rocs
  test <- lapply(tempList, 
                 function(x) roc.test(train_model_roc, x, method = comp_method)$p.value)
  # grab the actual AUC values
  auc_values <- t(as.data.frame.list(lapply(tempList, function(x) x$auc)))
  # Create a nice data frame to be outputed out
  aggregate_pvalues <- t(bind_cols(test)) %>% as.data.frame() %>% 
    mutate(study = rownames(.), BH = p.adjust(V1, method = "BH"), auc = auc_values[, 1]) %>% 
    rename(pvalue = V1) %>% select(study, auc, pvalue, BH)
  # Add the information on the test set
  aggregate_pvalues <- rbind(aggregate_pvalues, c(i, train_model_roc$auc, NA, NA))
  # return the final completed summary table
  return(aggregate_pvalues)
  
}


# Function to make comparisons between selected and full models
select_full_comparison <- function(full_model, select_model, 
                                   comp_method = "bootstrap"){
  # full_model is the model with all genera variables
  # select_model is the model with only crc specific variables
  # comp_method is default set to bootstrap to hedge against ROCs with different directions
  
  # Generates the pvalue from the test between the two respective models
  pvalue <- pROC::roc.test(full_model, select_model, 
                           method = comp_method)$p.value
  # creates a vector with auc or the two models and the pvalue
  all_data <- c(full_model = full_model$auc, select_model = select_model$auc, 
                pvalue = pvalue)
  # writes out the summary data
  return(all_data)
  
}

##############################################################################################
########################## Code used to run the analysis (unmatched) #########################
##############################################################################################

# reads in all the stool data into one list
unmatched_stool_study_data <- sapply(c(both_sets, tissue_sets), 
                                     function(x) get_data(x, tissue_unmatched), simplify = F)

#Align the genera so there is the same number for each data set.
unmatched_matched_genera_list <- align_genera(c(both_sets, tissue_sets), "column_length", 
                                    "sub_genera_data", unmatched_stool_study_data)

# Generate data sets to be used in random forest
unmatched_rf_datasets <- sapply(c(both_sets, tissue_sets), 
                      function(x) assign_disease(x, unmatched_matched_genera_list), simplify = F)

# Generate data for each test (study) set
unmatched_stool_final_data <- sapply(c(both_sets, tissue_sets), 
                           function(x) run_rf_tests(x, unmatched_rf_datasets), simplify = F)


# Generate summary data based on rocs
unmatched_pvalue_summaries <- sapply(names(unmatched_stool_final_data), 
                           function(x) make_model_comparisons(x, unmatched_stool_final_data), simplify = F)

# Generate final overal roc data for plotting
unmatched_all_roc_values <- sapply(names(unmatched_stool_final_data), 
                         function(x) make_data_table(unmatched_stool_final_data[[x]]), simplify = F)


##############################################################################################
########################## Code used to run the analysis (unmatched) #########################
##############################################################################################

# reads in all the stool data into one list
matched_stool_study_data <- sapply(c(tissue_sets), 
                                   function(x) get_data(x, tissue_matched), simplify = F)

#Align the genera so there is the same number for each data set.
matched_matched_genera_list <- list(lu = align_genera(c(tissue_sets), "column_length", 
                                              "sub_genera_data", matched_stool_study_data))

# Generate data sets to be used in random forest
matched_rf_datasets <- sapply(c(tissue_sets), 
                                function(x) assign_disease(x, matched_matched_genera_list), simplify = F)

# Generate data for each test (study) set
# Definitely overfit since the classification is 100%
matched_stool_model <- randomForest(disease ~ ., data = matched_rf_datasets[["lu"]], 
                                    mtry = round(sqrt(ncol(matched_rf_datasets[["lu"]]))), 
                                    importance = TRUE)

matched_model_rocs <- roc(matched_matched_genera_list$lu$disease ~ 
                            matched_stool_model$votes[, "polyp"])

# Get important OTUs to the model, are they relevant 
matched_model_importance_table <- matched_stool_model$importance %>% as.data.frame() %>% 
  mutate(genera = rownames(.)) %>% 
  arrange(desc(abs(MeanDecreaseAccuracy)))


##############################################################################################
############### Run the actual programs to get the data (CRC Specific Genera - unmatched) ####
##############################################################################################

# reduce the data sets down to only the CRC associated genera
select_unmatched_matched_genera_list <- lapply(unmatched_matched_genera_list, 
                                     function(x) 
                                       x %>% select(c("disease", "Fusobacterium", 
                                                      "Peptostreptococcus", 
                                                      "Porphyromonas", "Parvimonas")))

# Run the models
selected_unmatched_stool_final_data <- 
  sapply(c(both_sets, tissue_sets),  
  function(x) run_rf_tests(x, select_unmatched_matched_genera_list, specific_vars = T), simplify = F)


# Generate summary data based on rocs
selected_pvalue_summaries <- sapply(
  names(selected_unmatched_stool_final_data), 
  function(x) make_model_comparisons(x, selected_unmatched_stool_final_data), simplify = F)

# Generate final overal roc data for plotting
selected_unmatched_all_roc_values <- sapply(
  names(selected_unmatched_stool_final_data), 
  function(x) make_data_table(selected_unmatched_stool_final_data[[x]]), simplify = F)

# Compare the full data roc to the selected data roc and create a nice table
unmatched_test_red_select_models <- t(
  sapply(c(both_sets, tissue_sets), 
         function(x) 
           select_full_comparison(unmatched_stool_final_data[[x]][[x]], 
                                  selected_unmatched_stool_final_data[[x]][[x]]))) %>% 
  as.data.frame() %>% mutate(study = rownames(.), BH = p.adjust(pvalue, method = "BH")) %>% 
  select(study, full_model, select_model, pvalue, BH)

##############################################################################################
############### Run the actual programs to get the data (CRC Specific Genera - matched) ####
##############################################################################################

# reduce the data sets down to only the CRC associated genera
select_matched_matched_genera_list <- lapply(
  matched_matched_genera_list, 
  function(x) 
    x %>% 
    select(c("disease", "Fusobacterium", "Peptostreptococcus", "Porphyromonas", "Parvimonas")) %>% 
    mutate(disease = as.factor(disease)))


# Generate data for each test (study) set
# Definitely overfit since the classification is 100%
select_matched_stool_model <- randomForest(disease ~ ., data = select_matched_matched_genera_list[["lu"]], 
                                    mtry = round(sqrt(ncol(select_matched_matched_genera_list[["lu"]]))), 
                                    importance = TRUE)


select_matched_model_rocs <- roc(matched_matched_genera_list$lu$disease ~ 
                                   select_matched_stool_model$votes[, "polyp"])

# Compare the full data roc to the selected data roc and create a nice table
matched_test_red_select_models <- t(
  sapply(c(tissue_sets), 
         function(x) 
           select_full_comparison(matched_model_rocs, 
                                  select_matched_model_rocs))) %>% 
  as.data.frame() %>% mutate(study = rownames(.), BH = p.adjust(pvalue, method = "BH")) %>% 
  select(study, full_model, select_model, pvalue, BH)


















