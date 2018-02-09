### Run Random Forest Analysis -- Carcinoma
### Generate model and then test on remaining studies
### This code manually inserts data sets after the nzv step
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
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")


####################################################################################################
########################## Functions to run the analysis ###########################################
####################################################################################################


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
                         "stool", metadata = T) %>% filter(disease != "polyp")
  
  # Looks for Na in the meta data of interest and removes respective samples
  study_meta <- study_meta %>% filter(!is.na(disease))
  
  # conditional that checks for whether length of rows of meta data is smaller
  if(length(rownames(study_meta)) < length(rownames(sub_genera_data))){
    # grab only the samples in the meta data file for down stream analysis
    sub_genera_data <- sub_genera_data %>% slice(match(study_meta$sampleID, sample_ID))
    
    if(i == "hale"){
      
      study_meta <- study_meta %>% slice(match(sub_genera_data$sample_ID, sampleID))
    }
    
  } else{
    # grab only files in the data file for analysis
    study_meta <- study_meta %>% slice(match(sub_genera_data$sample_ID, sampleID))
  }
  # Prints out the total number of genera for that specific study
  print(paste("Total number of columns in", i, "is", 
              length(colnames(sub_genera_data))))
  # re assigns the rown names while removing the extra column used for sorting
  sample_names <- sub_genera_data$sample_ID
  # creates a list file with both data sets
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta, 
                   column_length = length(colnames(sub_genera_data)))
  # returns the combined list file
  return(dataList)
  
}


# Function that grabs the meta data and replaces sampleID with disease call
assign_disease <- function(studies, generaData, metaData, dataList){
  # studies is the variable with the names of the studies 
  # matched_genera is the data list with only genera in every study
  
  # Get the respective metadata file of interest
  tempData <- dataList[[studies]][[generaData]] %>% mutate(sample_ID = as.character(sample_ID))
  tempMeta <- dataList[[studies]][[metaData]] %>% mutate(sampleID = as.character(sampleID))
  # Gets transforms sample_ID column into a disease column with control v cancer calls
  corr_tempData <- tempData %>% left_join(select(tempMeta, sampleID, disease), by = c("sample_ID" = "sampleID")) %>% 
    mutate(disease = factor(ifelse(disease == "normal", 
                                   invisible("control"), invisible(disease)), 
                            levels = c("control", "cancer"))) %>% 
    select(-sample_ID) %>% 
    select(disease, everything())
  # Returns the modified data frame that can be used for RF analysis
  return(corr_tempData)
  
}



# Function to generate the training data sets based on removing near zero variance measures
create_training_data <- function(study, dataList){
  
  # gets the respective data set i for training
  training_data <- dataList[[study]]
  # stores the disease vector (it gets removed during processing for some studies)
  disease <- training_data$disease
  # Check for columns that have near zero variance
  nzv <- nearZeroVar(training_data)
  # remove columns that have near zero variance if nzv has value
  if(length(nzv) > 0){
    
    training_data <- training_data[, -nzv]
  }
  
  # Check to see if disease column was removed during the processing
  if("disease" %in% colnames(training_data)){
    # keep training data the same
    training_data <- training_data
    
  } else{
    # Re add disease to the training data at the beginning of the data table
    training_data <- training_data %>% mutate(disease = disease) %>% 
      select(disease, everything())
  }
  
  
  return(training_data)
}


# Function used to match genera and toss those in the test but not in the train
match_test_train_sets <- function(match_study, total_studies, dataList){
  
  train_genera <- colnames(dataList[[match_study]])
  
  tempList <- list()
  
  
  for(i in total_studies){
    
    if(match_study != i){
      
      temp_test_genera <- colnames(dataList[[i]])
      
      matched_genera <- c()
      unmatched_genera <- c()
      
      for(j in train_genera){
        
        if(j %in% temp_test_genera){
          
          matched_genera <- c(matched_genera, j)
        } else{
          
          unmatched_genera <- c(unmatched_genera, j)
        }
        
      }
      
      fake_data <- as.data.frame(sapply(unmatched_genera, 
                                        function(x) rep(0, length(rownames(dataList[[i]]))), 
                                        simplify = F) %>% bind_rows())
      
      tempList[[i]] <-  select(dataList[[i]], one_of(matched_genera)) %>% bind_cols(fake_data)
      
      
      
    }
  }
  
  return(tempList)
  
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
get_test_data <- function(train_study, i, 
                          training_model, training_data, testdataList){
  # i is the study of interest
  # train_study is the variable for the name of the study used as the training set 
  # training_model is the RF object of the training study
  # training_data is the raw data used to make the training model
  # testdataList is all the transformed data from all other studies
  
  # Get the outcomes of the voting for the training model
  train_prediction <- training_model[[train_study]]$finalModel$votes %>% as.data.frame()
  # get the predictions for each of the data sets based on the training model
  test_predictions <- sapply(i, function(x) 
    if(x != train_study) {
      predict(training_model[[train_study]], testdataList[[train_study]][[x]], 
              type = 'prob')}, simplify = F)
  
  # Generate roc curve infor (sens and spec) to be able to graph roc curves in the future
  overall_rocs <- sapply(i, function(x) 
    if(x != train_study){
      roc(testdataList[[train_study]][[x]]$disease ~ test_predictions[[x]][, "cancer"])}, 
    simplify = F)
  
  # add the training data roc information to this list
  overall_rocs[[train_study]] <- roc(training_data[[train_study]]$disease ~ train_prediction[, "cancer"])
  
  # Write out all the roc information from every data set
  return(overall_rocs)
  
}


# Function that will test all existing test sets (i.e. other studies)
get_select_test_data <- function(train_study, i, 
                                 training_model, testdataList){
  # i is the study of interest
  # train_study is the variable for the name of the study used as the training set 
  # training_model is the RF object of the training study
  # training_data is the raw data used to make the training model
  # testdataList is all the transformed data from all other studies
  
  # Get the outcomes of the voting for the training model
  train_prediction <- training_model[[train_study]]$finalModel$votes %>% as.data.frame()
  # get the predictions for each of the data sets based on the training model
  test_predictions <- sapply(i, function(x) 
    if(x != train_study) {
      predict(training_model[[train_study]], testdataList[[x]], 
              type = 'prob')}, simplify = F)
  
  # Generate roc curve infor (sens and spec) to be able to graph roc curves in the future
  overall_rocs <- sapply(i, function(x) 
    if(x != train_study){
      roc(testdataList[[x]]$disease ~ test_predictions[[x]][, "cancer"])}, 
    simplify = F)
  
  # add the training data roc information to this list
  overall_rocs[[train_study]] <- roc(testdataList[[train_study]]$disease ~ train_prediction[, "cancer"])
  
  # Write out all the roc information from every data set
  return(overall_rocs)
  
}


# Function to generate pvalues for the results ROC cuves (default method delong)
make_model_comparisons <- function(train_study, i, rocList, comp_method = "bootstrap"){
  # i is the study of interest
  # rocList is a list of roc objects (obtained from get_test_data)
  # comp_method is a character call of what method to use for comparisons
  # the default is set to bootstrap because it can test different direction curves
  
  # create a temp list variable 
  tempList <- rocList[[train_study]]
  # remove the study of interest (training roc) from temp list
  tempList[[train_study]][[train_study]] <- NULL
  # create the training roc variable
  train_model_roc <- rocList[[train_study]][[train_study]]
  # iterate through each study comparing the training roc the test rocs
  test <- lapply(tempList, function(x) 
    roc.test(train_model_roc, x, method = comp_method)$p.value)
  # grab the actual AUC values
  auc_values <- t(as.data.frame.list(lapply(tempList, function(x) x$auc)))
  # Create a nice data frame to be outputed out
  aggregate_pvalues <- t(bind_cols(test)) %>% as.data.frame() %>% 
    mutate(study = rownames(.), BH = p.adjust(V1, method = "BH"), auc = auc_values[, 1]) %>% 
    rename(pvalue = V1) %>% select(study, auc, pvalue, BH)
  # Add the information on the test set
  #aggregate_pvalues <- rbind(aggregate_pvalues, c(i, train_model_roc$auc, NA, NA))
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

# Function that gathers the important OTUs and takes the median with quartiles
get_imp_otu_data <- function(i, a_modelList){
  
  set.seed(12345)
  temp_rf_model <- make_rf_model(a_modelList[[i]])
  
  tempData <- varImp(temp_rf_model, scale = F)$importance %>% 
    as.data.frame() %>% 
    mutate(otu = rownames(.)) %>% 
    bind_rows() %>% 
    group_by(otu) %>% 
    summarise(mda_median = median(Overall), 
              iqr25 = quantile(Overall)["25%"], 
              iqr75 = quantile(Overall)["75%"]) %>% 
    arrange(desc(mda_median))
  
  return(tempData)
}



####################################################################################################
########################## Functions to run the analysis ###########################################
####################################################################################################

# reads in all the stool data into one list
stool_study_data <- mapply(get_data, c(stool_sets, "flemer"), SIMPLIFY = F)

# Generate data sets to be used in random forest
rf_datasets <- sapply(c(stool_sets, "flemer"), 
                      function(x) assign_disease(x, "sub_genera_data", "study_meta", stool_study_data), simplify = F)

training_imp_model_vars <- sapply(c(stool_sets, "flemer"), 
                                  function(x) get_imp_otu_data(x, rf_datasets), simplify = F)

# Generate training data
rf_training_data <- sapply(c(stool_sets, "flemer"), 
                           function(x) create_training_data(x, rf_datasets), simplify = F)

# Generate test data
rf_test_data <- sapply(c(stool_sets, "flemer"), 
                       function(x) match_test_train_sets(x, stool_sets, rf_training_data), simplify = F)


# Generate train models from training data
rf_training_models <- sapply(c(stool_sets, "flemer"), 
                             function(x) make_rf_model(rf_training_data[[x]]), simplify = F)


# Generate the data from testing on other studies
rf_study_test <- sapply(c(stool_sets, "flemer"), 
                        function(x) get_test_data(x, stool_sets, 
                                                  rf_training_models, rf_training_data, rf_test_data), simplify = F)

# Generate pvalue comparisons between train and test sets
train_test_pvalues <- sapply(c(stool_sets, "flemer"), 
                             function(x) make_model_comparisons(x, stool_sets, rf_study_test), simplify = F)



##############################################################################################
############### Run the actual programs to get the data (CRC Specific Genera) ################
##############################################################################################

rr_data <- read_csv("data/process/tables/select_genus_OR_stool_composite.csv") %>% arrange(pvalue, rr)

top5_pos_RR <- as.data.frame(rr_data %>% filter(rr > 1) %>% slice(1:5) %>% select(measure))[, "measure"]
top5_neg_RR <- as.data.frame(rr_data %>% filter(rr < 1) %>% slice(1:5) %>% select(measure))[, "measure"]
combined_genera <- c(top5_pos_RR, top5_neg_RR)

# reduce the data sets down to only the CRC associated genera
select_matched_genera_list <- sapply(names(stool_study_data), 
                                     function(x) 
                                       rf_datasets[[x]] %>% 
                                       select(c("disease", combined_genera)), simplify = F)

# Run the models
selected_rf_training_models <- sapply(
  c(stool_sets, "flemer"), 
  function(x) make_rf_model(select_matched_genera_list[[x]]), simplify = F)


# Generate the data from testing on other studies
selected_rf_study_test <- sapply(c(stool_sets, "flemer"), 
                                 function(x) get_select_test_data(x, stool_sets, 
                                                                  selected_rf_training_models, 
                                                                  select_matched_genera_list), simplify = F)


# Generate pvalue comparisons between train and test sets
selected_train_test_pvalues <- sapply(c(stool_sets, "flemer"), 
                                      function(x) make_model_comparisons(x, stool_sets, selected_rf_study_test), simplify = F)


# Compare the full data roc to the selected data roc and create a nice table
test_red_select_models <- sapply(stool_sets, 
                                 function(x) 
                                   as.data.frame(t(sapply(stool_sets, 
                                                          function(y) select_full_comparison(rf_study_test[[x]][[y]], 
                                                                                             selected_rf_study_test[[x]][[y]])))) %>% 
                                   mutate(study = rownames(.), train_model = x), 
                                 simplify = F) %>% bind_rows() %>% 
  mutate(BH = p.adjust(pvalue, method = "BH")) %>% 
  select(full_model, select_model, pvalue, BH, study, train_model)

##############################################################################################
############################## Write out the data ############################################
##############################################################################################

sapply(c(stool_sets, "flemer"), 
       function(x) write.csv(train_test_pvalues[[x]], 
                             paste("data/process/tables/ALL_genus_stool_RF_full_", 
                                   x, "_pvalue_summary.csv", sep = ""), row.names = F))

sapply(c(stool_sets, "flemer"), 
       function(x) write.csv(training_imp_model_vars[[x]], 
                             paste("data/process/tables/ALL_genus_stool_RF_full_", 
                                   x, "_imp_vars.csv", sep = ""), row.names = F))

sapply(c(stool_sets, "flemer"), 
       function(x) write.csv(selected_train_test_pvalues[[x]], 
                             paste("data/process/tables/ALL_genus_stool_RF_select_", 
                                   x, "_pvalue_summary.csv", sep = ""), row.names = F))

write.csv(test_red_select_models, 
          "data/process/tables/ALL_genus_stool_RF_fullvsselect_pvalue_summary.csv", row.names = F)














