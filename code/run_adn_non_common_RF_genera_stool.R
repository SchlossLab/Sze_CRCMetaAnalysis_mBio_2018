### Run Random Forest Analysis -- Adenoma
### Generate model and then test on remaining studies
### This code manually inserts data sets after the nzv step
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "caret", "pROC"))


# Stool Only polyp sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("brim", "zeller", "baxter", "hale")


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
    select(sample_ID, everything()) %>% mutate(sample_ID = as.character(sample_ID))
  # grabs the meta data and transforms polyp to control (polyp/control vs cancer) 
  study_meta <- get_file(i, "data/process/", ".metadata", rows_present = F,  
                         "stool", metadata = T) %>% 
    filter(disease != "cancer", !is.na(disease)) %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    select(sampleID, disease)
  
  sub_genera_data <- study_meta %>% 
    inner_join(sub_genera_data, by = c("sampleID" = "sample_ID")) %>% 
    select(-sampleID)
  
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta, 
                   column_length = length(colnames(sub_genera_data)))
  # returns the combined list file
  return(dataList)
  
}


# Function that grabs the meta data and replaces sampleID with disease call
assign_disease <- function(studies, generaData, dataList){
  # studies is the variable with the names of the studies 
  # matched_genera is the data list with only genera in every study
  
  # Get the respective metadata file of interest
  tempData <- dataList[[studies]][[generaData]]
  # Gets transforms sample_ID column into a disease column with control v cancer calls
  corr_tempData <- tempData %>% 
    mutate(disease = factor(ifelse(disease == "normal", 
                                   invisible("control"), invisible(disease)), 
                            levels = c("control", "polyp")))
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
      roc(testdataList[[train_study]][[x]]$disease ~ test_predictions[[x]][, "polyp"])}, 
    simplify = F)
  
  # add the training data roc information to this list
  overall_rocs[[train_study]] <- roc(training_data[[train_study]]$disease ~ train_prediction[, "polyp"])
  
  # Write out all the roc information from every data set
  return(overall_rocs)
  
}







####################################################################################################
########################## Functions to run the analysis ###########################################
####################################################################################################

# reads in all the stool data into one list
stool_study_data <- mapply(get_data, stool_sets, SIMPLIFY = F)

# Generate data sets to be used in random forest
rf_datasets <- sapply(stool_sets, 
                      function(x) assign_disease(x, "sub_genera_data", stool_study_data), simplify = F)

# Generate the test sets
rf_test_sets <- sapply(stool_sets, 
                       function(x) create_training_data(x, rf_datasets), simplify = F)

# Generate training data
rf_training_data <- sapply(stool_sets, 
                           function(x) create_training_data(x, rf_datasets), simplify = F)

# Generate test data
rf_test_data <- sapply(stool_sets, 
                       function(x) match_test_train_sets(x, stool_sets, rf_training_data), simplify = F)


# Generate train models from training data
rf_training_models <- sapply(stool_sets, 
                             function(x) make_rf_model(rf_training_data[[x]]), simplify = F)


# Generate the data from testing on other studies
rf_study_test <- sapply(stool_sets, 
                        function(x) get_test_data(x, stool_sets, 
                      rf_training_models, rf_training_data, rf_test_data), simplify = F)




