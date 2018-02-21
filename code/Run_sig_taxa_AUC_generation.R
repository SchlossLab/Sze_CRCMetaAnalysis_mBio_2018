### Generate counts and AUCs from significant OR data 
### Stool
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "caret", "pROC"))

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


# Function to generate step-wise increase in thresholds
get_thresholds <- function(study, col_of_int, dataList){
  
  tempData <- dataList[[study]] %>% 
    select(one_of(col_of_int))
  
  tempMax <- apply(tempData, 2, function(x) max(x))
  
  
  return(tempMax)
  
}


test <- get_thresholds("baxter", combined_genera, select_good_data)

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


