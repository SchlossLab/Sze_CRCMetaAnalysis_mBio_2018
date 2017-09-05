### Pull, transform, and normalize 4 crc genera 
### Specifically analyze the genera tied to crc from previous research
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "vegan", "foreach", "doParallel"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
tissue_sets <- c("dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "flemer", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")

# CRC genera
crc_genera <- c("Fusobacterium", "Peptostreptococcus", "Porphyromonas", "Parvimonas")



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
  # re assigns the rown names while removing the extra column used for sorting
  sample_names <- sub_genera_data$sample_ID
  # creates a list file with both data sets
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta)
  # returns the combined list file
  return(dataList)
  
}


# Function to grab only the genera file and pull specific genus from it
get_specific_genera <- function(i, genera_to_get, table_name, meta_name, 
                                dataList = stool_study_data){
  
  tempData <- dataList[[i]][[table_name]] %>% 
    select(sample_ID, one_of(genera_to_get)) %>% 
    rename(sampleID = sample_ID) %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    inner_join((dataList[[i]][[meta_name]] %>% 
                  mutate(sampleID = as.character(sampleID), 
                         disease = ifelse(disease == "normal", 
                                          invisible("control"), invisible(disease)))), 
               by = "sampleID") %>% 
    mutate(disease2 = ifelse(disease != "cancer", invisible("control"), invisible(disease))) %>% 
      as.data.frame()
  
  print(paste("Finished pulling and merging files for", i, "data sets"))
  
  return(tempData)
  
}



# Function to judge whether lambda is negative or positive and obtain proper lambda
calc_lambda_sign <- function(filtered_vec, transformed_vec){
  # filtered_vec are original alpha diversity values
  # transformed_vec are the appropriate power transformed values
  
  # Runs to check if the transformed value is positive or negative
  if(transformed_vec > 0){
    # Conversion formula if value is positive
    lambda_value <- log(transformed_vec)/log(filtered_vec)
  } else{
    # Conversion formula if value is negative
    lambda_value <- log(-1*transformed_vec)/log(filtered_vec)
  }
  # Write out the obtained value
  return(lambda_value)
}


# Control function to run through each alpha metric column of supplied data frames
generate_lambda <- function(i, original_df, transformed_data){
  # i is a vector of alpha metrics
  # original_df is a data frame with both meta data and alpha metrics
  # transformed_data is a data frame with power transformed alpha metrics of interest
  
  # create filtered data frame with only alpha metrics of interest
  filtered_df <- original_df %>% 
    select(sobs, shannon, shannoneven)
  
  # converts provided data frames to a data frame (tibbles make things screwy)
  conv_filtered_data <- as.data.frame(filtered_df[1, i])
  conv_transformed_data <- as.data.frame(transformed_data[1, i])
  
  # run the calculation function on every alpha metric of interest
  calculated_lambda <- mapply(calc_lambda_sign, conv_filtered_data, conv_transformed_data)
  
  # Write out the obtained values
  return(unname(calculated_lambda))
}


# Function to transform needed data, get the lambda and p-values for normality
get_lambda <- function(test_df){
  # test_df is a data frame with metadata and alpha metrics of interest
  
  # Apply optimal power transformation
  new_power_data <- test_df %>% 
    mutate_at(vars(sobs, shannon, shannoneven), 
              function(x) transformTukey(x, plotit = FALSE))
  
  test_new_data <- new_power_data %>% select(sobs, shannon, shannoneven)
  
  optimum_lambdas <- generate_lambda(c("sobs", "shannon", "shannoneven"), 
                                     test_df, test_new_data)
  
  pvalues <- apply(test_new_data, 2, function(x) shapiro.test(x)$p.value)
  
  values <- cbind(optimum_lambdas, pvalues)
  
  lambda_data <- list(summary_stats = values, transformed_data = new_power_data)
  
  return(lambda_data)
  
}






##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# reads in all the stool data into one list
stool_study_data <- mapply(get_data, stool_sets, SIMPLIFY = F)


specific_genera_list <- sapply(c(stool_sets, "flemer"), 
               function(x) get_specific_genera(x, crc_genera, 
                                               "sub_genera_data", "study_meta"))











