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
  
  return(tempData)
  
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


train_data <- rf_datasets[["wang"]]

nzv <- nearZeroVar(train_data)

train_data <- train_data[, -nzv]

tempList <- rf_datasets

tempList[["wang"]] <- NULL

tempList <- lapply(tempList, function(x) as.data.frame(x[, -nzv]))


preProcValues <- preProcess(train_data, method = c("YeoJohnson", "center", "scale"))

train_data_transformed <- predict(preProcValues, train_data)

tempList <- lapply(tempList, function(x) predict(preProcValues, x))

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
train_model_data <- 
  train(disease ~ ., data = train_data, 
        method = "rf", 
        ntree = 500, 
        trControl = fitControl,
        tuneGrid = tunegrid, 
        metric = "ROC", 
        na.action = na.omit, 
        verbose = FALSE)

# Generate data for each test (study) set
test <- get_test_data(names(tempList), "wang", train_model_data, train_data, tempList)






# Set up the pulling of important information

test_sens <- test_roc$sensitivities
test_spec <- test_roc$specificities
test_auc <- ifelse(test_roc$auc < 0.5, 
                    invisible(1-test_roc$auc), invisible(test_roc$auc))


test2 <- sapply(names(test), function(x) 
  as.data.frame(cbind(sens = as.numeric(test[[x]]$sensitivities), 
        spec = test[[x]]$specificities, 
        auc = rep(ifelse(
          test[[x]]$auc[1] < 0.5, invisible(1 - test[[x]]$auc[1]), 
          invisible(test[[x]]$auc[1])), 
          length(test[[x]]$sensitivities)), 
        study = rep(x, length(test[[x]]$sensitivities))), stringsAsFactors = F), simplify = F)


train_sens <- train_roc$sensitivities
train_spec <- train_roc$specificities
train_mtry <- train_model_data$results$mtry
train_auc <- ifelse(train_model_data$results$ROC < 0.5, 
                    invisible(1-train_model_data$results$ROC), 
                    invisible(train_model_data$results$ROC))











#### TO DO LIST ####
#### Set up a function to iterate and store all possible combinations

