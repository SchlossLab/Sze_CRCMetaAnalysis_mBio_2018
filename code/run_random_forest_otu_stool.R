### Run Random Forest Analysis OTU Level
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
  
  # grabs subsampled data and assigns rownames from sample names to table
  shared_data <- read.delim(paste("data/process/", i, "/", i, ".0.03.subsample.shared", 
                                  sep = ""), header = T, stringsAsFactors = F) %>% 
    select(-label, -numOtus)

  # grabs the meta data and transforms polyp to control (polyp/control vs cancer) 
  study_meta <- get_file(i, "data/process/", ".metadata", rows_present = F,  
                         "stool", metadata = T)
  
  # Looks for Na in the meta data of interest and removes respective samples
  study_meta <- study_meta %>% filter(!is.na(disease))
  
  # conditional that checks for whether length of rows of meta data is smaller
  if(length(rownames(study_meta)) < length(rownames(shared_data))){
    # grab only the samples in the meta data file for down stream analysis
    shared_data <- shared_data %>% slice(match(study_meta$sampleID, Group))
    
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
  
  if(length(nzv) == 0){
    
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





## TO Do List

## Run a 10-fold CV (5 if 10 is not possible)

## Run a random label 10-fold CV (5 if 10 is not possible)

## Compare the two models

## Save the ROC curves to graph 

## Introduce For each loop possibly






##############################################################################################
############### Run the actual programs to get the data (ALL Data) ###########################
##############################################################################################


for(i in "weir"){
  
  dataList <- get_data(i = i)
  
  disease_dataset <- assign_disease("study_meta", "shared_data", dataList)
  
  rf_data <- get_align_info(disease_dataset)
  
  actual_model <- make_rf_model(rf_data[["train_data"]])
  
  random_model <- make_rf_model(rf_data[["rand_data"]])
}
















