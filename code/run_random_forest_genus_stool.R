### Run Random Forest Analysis
### Generate model and then test on remaining studies
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "caret", "pROC"))

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
                         "stool", metadata = T)
  
  # Looks for Na in the meta data of interest and removes respective samples
  study_meta <- study_meta %>% filter(!is.na(disease))
  
  # conditional that checks for whether length of rows of meta data is smaller
  if(length(rownames(study_meta)) < length(rownames(sub_genera_data))){
    # grab only the samples in the meta data file for down stream analysis
    sub_genera_data <- sub_genera_data %>% slice(match(study_meta$sampleID, sample_ID))
    
  } else{
    # grab only files in the data file for analysis
    study_meta <- study_meta %>% slice(match(sub_genera_data$sample_ID, sampleID))
  }
  
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
  
  genera_num_list <- sort.int(sapply(studies, 
                                     function(x) dataList[[x]][[length_column_name]]))
  
  x = 1
  
  while(genera_num_list[1]!= genera_num_list[length(genera_num_list)]){
    
    lowest_genera_study <- names(genera_num_list[1])
    
    if(x == 1){
      
      genera_names <- dataList[[lowest_genera_study]][[genera_data_name]] %>% 
        select(-contains("_unclassified")) %>% colnames(.)
    } else{
      
      genera_names <- colnames(temp_aligned_genera[[lowest_genera_study]])
    }
    
    
    temp_aligned_genera <- suppressWarnings(sapply(studies, 
                                  function(x) 
                                    select(dataList[[x]][[genera_data_name]], 
                                           one_of(genera_names)), simplify = F))
    
    genera_num_list <- sort.int(sapply(studies, 
                                    function(x) 
                                      length(colnames(temp_aligned_genera[[x]]))))
    
    print(paste("Min and Max total genera is:", 
                min(genera_num_list), ",", max(genera_num_list)))
    
    x = x + 1
    
  }

  return(temp_aligned_genera)
}


# Function that grabs the meta data and replaces sampleID with disease call
assign_disease <- function(studies, metadata_table_name, 
                           matched_genera, fullDataList){
  
  tempMetadata <- fullDataList[[studies]][[metadata_table_name]]
  
  tempData <- matched_genera[[studies]] %>% 
    mutate(sample_ID = factor(ifelse(tempMetadata$disease == "polyp", 
                                     invisible("control"), 
                                     ifelse(tempMetadata$disease == "normal", 
                                            invisible("control"), 
                                            invisible(tempMetadata$disease))), 
                              levels = c("control", "cancer"))) %>% 
    rename(disease = sample_ID)
  
  return(as.data.frame(tempData))
  
}


# Function to apply and get the nzv and preProcess for the training data
get_align_info <- function(i, dataList){
  
  training_data <- dataList[[i]]
  
  disease <- training_data$disease
  
  nzv <- nearZeroVar(training_data)
  
  training_data <- training_data[, -nzv]
  
  preProcValues <- preProcess(training_data, method = c("YeoJohnson", "center", "scale"))
  
  training_data <- predict(preProcValues, training_data)
  
  if("disease" %in% colnames(training_data)){
    
    training_data <- training_data
    
  } else{
    
    training_data <- training_data %>% mutate(disease = disease) %>% 
      select(disease, everything())
  }
  
  final_info <- list(train_data = training_data, 
                     near_zero_variance = nzv, 
                     scaling = preProcValues)

  return(final_info)
}


# Function to generate the preprocessing files on test data
apply_preprocess <- function(i, trainingList_info, dataList){
  
  dataList[[i]] <- NULL
  
  nzv <- as.numeric(trainingList_info[["near_zero_variance"]])
  
  if(1 %in% nzv){
    
    nzv <- nzv[nzv != 1]
    
  } else {
    
    nzv <- nzv
  }
  
  preProcValues <- trainingList_info[["scaling"]]
  
  dataList <- lapply(dataList, function(x) as.data.frame(x[, -nzv]))
  dataList <- lapply(dataList, function(x) predict(preProcValues, x))
  
  return(dataList)
  
}


# Function that will run and create the needed model
make_rf_model <- function(train_data){
  
  #Create Overall specifications for model tuning
  # number controls fold of cross validation
  # Repeats control the number of times to run it
  
  fitControl <- trainControl(## 5-fold CV
    method = "cv",
    number = 10,
    p = 0.8, 
    classProbs = TRUE, 
    summaryFunction = twoClassSummary, 
    savePredictions = "final")
  
  number_try <- round(sqrt(ncol(train_data)))
  
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
  
  return(training_model)
}


# Function that will test all existing test sets (i.e. other studies)
get_test_data <- function(i, train_study, 
                          training_model, training_data, testdataList){
  
  # Set up testing
  train_prediction <- training_model$finalModel$votes %>% as.data.frame()
  
  test_predictions <- sapply(i, function(x) 
    predict(training_model, testdataList[[x]], type = 'prob'), simplify = F)
    
  overall_rocs <- sapply(i, function(x) 
    roc(testdataList[[x]]$disease ~ test_predictions[[x]][, "cancer"]), simplify = F)
  
  overall_rocs[[train_study]] <- roc(training_data$disease ~ train_prediction[, "cancer"])
  

  return(overall_rocs)
  
}

# Function that creates a finalized data table for respective model
make_data_table <- function(final_rocs){
  
  tempData <- sapply(names(final_rocs), function(x) 
    as.data.frame(cbind(sens = as.numeric(final_rocs[[x]]$sensitivities), 
                        spec = final_rocs[[x]]$specificities, 
                        auc = rep(final_rocs[[x]]$auc[1], length(final_rocs[[x]]$sensitivities))), 
                  stringsAsFactors = F) %>% 
      mutate(study = rep(x, length(final_rocs[[x]]$sensitivities))), simplify = F) %>% bind_rows()
  
  return(tempData)
}


# Function to execute the major commands for RF gathering
run_rf_tests <- function(study, rf_dataList){
  
  # Generate data for each test (study) set
  first_study <- get_align_info(study, rf_dataList)
  
  test_dataList <- apply_preprocess(study, first_study, rf_dataList)
  
  
  train_model_data <- make_rf_model(first_study$train_data)
  
  
  test_dataLists <- get_test_data(names(test_dataList), study, train_model_data, 
                        first_study[["train_data"]], test_dataList)
  
  final_results <- make_data_table(test_dataLists)
  
  return(final_results)
}



##############################################################################################
############### Run the actual programs to get the data ######################################
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
final_data <- sapply(c(stool_sets, "flemer"), 
                     function(x) run_rf_tests(x, rf_datasets), simplify = F)












#### TO DO LIST ####
#### Set up a function to iterate and store all possible combinations

