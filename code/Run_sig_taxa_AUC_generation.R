### Generate counts and AUCs from significant OR data 
### Stool
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret", "pROC", "AUC"))

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "hale")

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


# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_tissue_data <- function(i, metadata){
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
  study_meta <- metadata %>% filter(study == i)
  
  # Looks for Na in the meta data of interest and removes respective samples
  study_meta <- study_meta %>% filter(!is.na(disease))
  
  # conditional that checks for whether length of rows of meta data is smaller
  if(length(rownames(study_meta)) < length(rownames(sub_genera_data))){
    # grab only the samples in the meta data file for down stream analysis
    sub_genera_data <- sub_genera_data %>% slice(match(study_meta$sample_id, sample_ID))
    
  } else{
    # grab only files in the data file for analysis
    study_meta <- study_meta %>% slice(match(sub_genera_data$sample_ID, sample_id))
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
assign_disease <- function(studies, generaData, metaData, dataList, stool = T){
  # studies is the variable with the names of the studies 
  # matched_genera is the data list with only genera in every study
  
  if(stool == T){
    
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
    
  } else{
    # Get the respective metadata file of interest
    tempData <- dataList[[studies]][[generaData]] %>% mutate(sample_ID = as.character(sample_ID))
    tempMeta <- dataList[[studies]][[metaData]] %>% mutate(sample_ID = as.character(sample_id))
    # Gets transforms sample_ID column into a disease column with control v cancer calls
    corr_tempData <- tempData %>% 
      left_join(select(tempMeta, sample_id, disease), by = c("sample_ID" = "sample_id")) %>% 
      filter(disease != "polyp") %>% 
      mutate(disease = factor(ifelse(disease == "normal", 
                                     invisible("control"), invisible(disease)), 
                              levels = c("control", "cancer"))) %>% 
      select(-sample_ID) %>% 
      select(disease, everything())
  }
  
  
  # Returns the modified data frame that can be used for RF analysis
  return(corr_tempData)
  
}


# Function to generate step-wise increase in thresholds
get_thresholds <- function(study, col_of_int, dataList){
  
  tempData <- dataList[[study]] %>% 
    select(one_of(col_of_int))
  
  actualIntervals <- apply(tempData, 2, function(x) median(x))
  
  return(actualIntervals)
  
}

# Function to generate counts (TP, TN, FP, FN) from thresholds
generate_counts <- function(study, col_of_int, thresholdList, dataList){
  
  tempData <- as.data.frame(dataList[[study]])
  tempThreshold <- thresholdList[[study]]
  tempDisease <- tempData$disease
  
  
  tempGeneraAnalysis <- sapply(col_of_int, 
                               function(x) making_counts(x, tempThreshold, tempData), simplify = F)
  
  
  test_counts <- t(as.data.frame.list(sapply(col_of_int, 
                        function(x) get_tp_tn(x, tempDisease, tempGeneraAnalysis)))) %>% 
    as.data.frame() %>% 
    mutate(temp_col = rownames(.),
           study = study) %>% 
    separate(temp_col, c("taxa", "e1, e2"), sep = "[.]") %>% 
    rename(auc = V1) %>% 
    select(taxa, auc, study)
  
  return(test_counts)
  
}


# Function to generate the disease calls from thresholds
making_counts <- function(taxa, thresholds, dataTable){
  
  tempVector <- dataTable[, taxa]
  tempThreshold <- thresholds[taxa]
  
  tempTest <- sapply(1:length(tempThreshold), 
                     function(x) ifelse(tempVector < tempThreshold[x], 
                                        invisible("control"), invisible("cancer")), simplify = F)
  
  
  return(tempTest)
  
}

# Function to generate the actual counts
get_tp_tn <- function(taxa, diseaseVector, dataList){
  
  tempList <- dataList[[taxa]]
  tempDataStore <- list()
  
  tempAUC <- sapply(1:length(tempList), 
                    function(x) auc(roc(tempList[[x]], diseaseVector)), simplify = F)
  

  return(tempAUC)
  
}



##############################################################################################
############### Run the actual programs to get the data (CRC Specific Genera) ################
##############################################################################################

# reads in all the stool data into one list
stool_study_data <- mapply(get_data, c(stool_sets, "flemer"), SIMPLIFY = F)

# Generate data sets to be used in random forest
good_datasets <- sapply(c(stool_sets, "flemer"), 
                      function(x) assign_disease(x, "sub_genera_data", "study_meta", stool_study_data), simplify = F)

# Grab the significant taxa
rr_data <- read_csv("data/process/tables/select_genus_OR_stool_composite.csv") %>% 
  arrange(pvalue, rr) %>% 
  mutate(bh = p.adjust(pvalue, method = "BH")) %>% 
  filter(bh < 0.05) %>% 
  select(measure)

combined_genera <- rr_data$measure

# pick out only significant taxa from each data set
select_good_data <- sapply(c(stool_sets, "flemer"), 
                           function(x) select(good_datasets[[x]], 
                                              disease, one_of(combined_genera)), simplify = F)

# Get the thresholds for each study
study_thresholds <- sapply(c(stool_sets, "flemer"), 
                           function(x) get_thresholds(x, combined_genera, select_good_data), simplify = F)


# Generate the counts for each study and each taxa
count_data_list <- sapply(c(stool_sets, "flemer"), 
                          function(x) generate_counts(x, combined_genera, study_thresholds, select_good_data), 
                          simplify = F) %>% bind_rows()


# Read out final data
write_csv(count_data_list, "data/process/tables/ind_genera_auc_stool.csv")


##############################################################################################
############### Run the actual programs to get the data (CRC Unmatched Tissue) ################
##############################################################################################

# Read in specific data tables to be used
tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  mutate(disease = gsub("adenoma", "polyp", disease)) %>% 
  select(group, id, disease, study) %>% 
  rename(sample_id = group)

# Get studies that contain unmatched samples
# Need to remove dejea and lu
# Lu only has polyps
# Dejea only has cancer 
unmatched_studies <- unique(
  tissue_unmatched$study[!(tissue_unmatched$study %in% c("dejea", "lu"))]) 

# reads in all the stool data into one list
unmatched_study_data <- sapply(unmatched_studies, 
                                     function(x) get_tissue_data(x, tissue_unmatched), simplify = F)

# Generate data sets to be used in random forest
unmatched_good_datasets <- sapply(unmatched_studies, 
                        function(x) assign_disease(x, "sub_genera_data", "study_meta", unmatched_study_data, stool = F), 
                        simplify = F)

rr_data <- read_csv("data/process/tables/select_genus_OR_unmatched_tissue_composite.csv") %>% 
  arrange(pvalue, rr) %>% 
  mutate(bh = p.adjust(pvalue, method = "BH")) %>% 
  filter(bh < 0.05) %>% 
  select(measure)

combined_genera <- rr_data$measure

# pick out only significant taxa from each data set
select_tissue_good_data <- sapply(unmatched_studies, 
                           function(x) select(unmatched_good_datasets[[x]], 
                                              disease, one_of(combined_genera)), simplify = F)

# Get the thresholds for each study
study_tissue_thresholds <- sapply(unmatched_studies, 
                           function(x) get_thresholds(x, combined_genera, select_tissue_good_data), simplify = F)


# Generate the counts for each study and each taxa
count_tissue_data_list <- sapply(unmatched_studies, 
                          function(x) generate_counts(x, combined_genera, study_tissue_thresholds, select_tissue_good_data), 
                          simplify = F) %>% bind_rows()


# Read out final data
write_csv(count_data_list, "data/process/tables/ind_genera_auc_unmatched_tissue.csv")

