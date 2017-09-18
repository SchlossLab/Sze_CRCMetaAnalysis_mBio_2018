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
  
  
}
















