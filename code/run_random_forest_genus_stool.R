### Run Random Forest Analysis
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
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")


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


# Function that aligns the genera from all data sets
align_genera <- function(studies, length_column_name, 
                         genera_data_name, dataList){
  # studies is a variable for the name of study to be used
  # length_column_name is the name of the vector in dataList that stores the total genera number
  # genera_data_name is the name of the data table in dataList that contains the genera information
  # dataList is a list that has for every study the genus info, metadata, and total genera present
  
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
  # When the while loop exits return the aligned genera data list
  return(temp_aligned_genera)
}


# Function that grabs the meta data and replaces sampleID with disease call
assign_disease <- function(studies, metadata_table_name, 
                           matched_genera, fullDataList){
  # studies is the variable with the names of the studies 
  # metadata_table_name is the variable with the name of th metadata file
  # matched_genera is the data list with only genera in every study
  # fullDataList is the original created data list
  
  # Get the respective metadata file of interest
  tempMetadata <- fullDataList[[studies]][[metadata_table_name]]
  # Gets transforms sample_ID column into a disease column with control v cancer calls
  tempData <- matched_genera[[studies]] %>% 
    mutate(sample_ID = factor(ifelse(tempMetadata$disease == "polyp", 
                                     invisible("control"), 
                                     ifelse(tempMetadata$disease == "normal", 
                                            invisible("control"), 
                                            invisible(tempMetadata$disease))), 
                              levels = c("control", "cancer"))) %>% 
    rename(disease = sample_ID)
  # Returns the modified data frame that can be used for RF analysis
  return(as.data.frame(tempData))
  
}


# Function to apply and get the nzv and preProcess for the training data
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
  # Find the Power transform, center, and zscore normalize formula
  preProcValues <- preProcess(training_data, method = c("YeoJohnson", "center", "scale"))
  # Perform the transformation on the training data
  training_data <- predict(preProcValues, training_data)
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
                     near_zero_variance = nzv, 
                     scaling = preProcValues)
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
  # pull out the formula/object to be used for transforming all the data
  preProcValues <- trainingList_info[["scaling"]]
  # Remove near zero variance from all data sets
  dataList <- lapply(dataList, function(x) as.data.frame(x[, -nzv]))
  # Transform all the data using the same transformations applied to the training data
  dataList <- lapply(dataList, function(x) predict(preProcValues, x))
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
    roc(testdataList[[x]]$disease ~ test_predictions[[x]][, "cancer"]), simplify = F)
  # add the training data roc information to this list
  overall_rocs[[train_study]] <- roc(training_data$disease ~ train_prediction[, "cancer"])
  
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
  
  if(specific_vars == F){
    
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


# Function that gathers the important OTUs and takes the median with quartiles
get_imp_otu_data <- function(i, a_modelList){
  
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



##############################################################################################
############### Run the actual programs to get the data (ALL Data) ###########################
##############################################################################################

# reads in all the stool data into one list
stool_study_data <- mapply(get_data, c(stool_sets, "flemer"), SIMPLIFY = F)

#Align the genera so there is the same number for each data set.
matched_genera_list <- align_genera(c(stool_sets, "flemer"), "column_length", 
                                    "sub_genera_data", stool_study_data)

# Generate data sets to be used in random forest
rf_datasets <- sapply(c(stool_sets, "flemer"), 
                      function(x) assign_disease(x, "study_meta", 
                                                 matched_genera_list, stool_study_data), 
                      simplify = F)


# Generate data for each test (study) set
stool_final_data <- sapply(c(stool_sets, "flemer"), 
                     function(x) run_rf_tests(x, rf_datasets), simplify = F)


# generate the important variables in the model
imp_vars <- sapply(c(stool_sets, "flemer"), 
                   function(x) get_imp_otu_data(x, rf_datasets), simplify = F)

# Generate summary data based on rocs
pvalue_summaries <- sapply(names(stool_final_data), 
               function(x) make_model_comparisons(x, stool_final_data), simplify = F)

# Generate final overal roc data for plotting
all_roc_values <- sapply(names(stool_final_data), 
                   function(x) make_data_table(stool_final_data[[x]]), simplify = F)



##############################################################################################
############### Run the actual programs to get the data (CRC Specific Genera) ################
##############################################################################################

# reduce the data sets down to only the CRC associated genera
select_matched_genera_list <- lapply(matched_genera_list, 
                                     function(x) 
                                       x %>% select(c("sample_ID", "Fusobacterium", 
                                                      "Peptostreptococcus", 
                                                      "Porphyromonas", "Parvimonas")))

# Generate data sets to be used in random forest
selected_rf_datasets <- sapply(c(stool_sets, "flemer"), 
                      function(x) 
                        assign_disease(x, "study_meta", 
                                       select_matched_genera_list, stool_study_data), simplify = F)

# Run the models
selected_stool_final_data <- sapply(c(stool_sets, "flemer"), 
                           function(x) 
                             run_rf_tests(x, selected_rf_datasets, specific_vars = T), simplify = F)


# generate the important variables in the model
selected_imp_vars <- sapply(c(stool_sets, "flemer"), 
                   function(x) get_imp_otu_data(x, rf_datasets), simplify = F)

# Generate summary data based on rocs
selected_pvalue_summaries <- sapply(names(selected_stool_final_data), 
                           function(x) 
                             make_model_comparisons(x, selected_stool_final_data), simplify = F)

# Generate final overal roc data for plotting
selected_all_roc_values <- sapply(names(selected_stool_final_data), 
                         function(x) make_data_table(selected_stool_final_data[[x]]), simplify = F)

# Compare the full data roc to the selected data roc and create a nice table
test_red_select_models <- t(sapply(c(stool_sets, "flemer"), 
               function(x) 
                 select_full_comparison(stool_final_data[[x]][[x]], 
                                        selected_stool_final_data[[x]][[x]]))) %>% 
  as.data.frame() %>% mutate(study = rownames(.), BH = p.adjust(pvalue, method = "BH")) %>% 
  select(study, full_model, select_model, pvalue, BH)


##############################################################################################
############################## Write out the data ############################################
##############################################################################################

# Write out all the necessary data table files

sapply(c(stool_sets, "flemer"), 
       function(x) write.csv(pvalue_summaries[[x]], 
                             paste("data/process/tables/genus_stool_RF_full_", 
                                   x, "_pvalue_summary.csv", sep = ""), row.names = F))

sapply(c(stool_sets, "flemer"), 
       function(x) write.csv(all_roc_values[[x]], 
                             paste("data/process/tables/genus_stool_RF_full_", 
                                   x, "_raw_roc_data.csv", sep = ""), row.names = F))

sapply(c(stool_sets, "flemer"), 
       function(x) write.csv(selected_pvalue_summaries[[x]], 
                             paste("data/process/tables/genus_stool_RF_select_", 
                                   x, "_pvalue_summary.csv", sep = ""), row.names = F))


sapply(c(stool_sets, "flemer"), 
       function(x) write.csv(selected_all_roc_values[[x]], 
                             paste("data/process/tables/genus_stool_RF_select_", 
                                   x, "_raw_roc_data.csv", sep = ""), row.names = F))

write.csv(test_red_select_models, 
          "data/process/tables/genus_stool_RF_fullvsselect_pvalue_summary.csv", row.names = F)














