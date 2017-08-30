### Get the dmm groups 
### Merge the group data with relevant metadata
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
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "flemer")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")

# The package DirichletMultinomial really messes with dplyr and makes it almost impossible to
# use things such as slice and match in combination. Make sure to do all table 
# transformations before loading and using it.


# Function to read in needed genera file
get_file <- function(i, path_to_file, ending, rows_present=T, name_of_rows=1, 
                     vec_of_rownames = NULL, sample_source, metadata = F){
  # i represents the study name
  # path_to_file is the is the directory where the main studies are stored
  # ending is the extension after i of the name of the file
  # rows_present signifies whether rownames are part of the inputed file
  # name_of_rows represents what column the row names are in
  # sample_source signifies whether it is stool or tissue
  # metadata signifies whether the inputed file is metadata or not
  
  # First conditional to make sure only data goes through this part
  if(metadata == F){
    # second conditional to check if rownames are present or not
    if(rows_present == T){
      # reads in data (assumes csv file) with row names from data table
      temp_shared <- read.csv(paste(path_to_file, i, "/", i, ending, sep = ""), 
                              header = T, stringsAsFactors = F, row.names = name_of_rows)
    } else{
      # reads in data (assumes csv file) with row names from an inputed vector 
      temp_shared <- read.csv(paste(path_to_file, i, "/", i, ending, sep = ""), 
                              header = T, stringsAsFactors = F)
      rownames(temp_shared) <- vec_of_rownames
    }
    
  } else{
    # reads in data table (assumes tab delimited) for meta data
    temp_shared <- read.delim(paste(path_to_file, i, "/", i, ending, sep = ""), 
                              header = T, stringsAsFactors = F)
    
    # conditional that checks for column named sample_type (created during mothur processing)
    if(!("sample_type" %in% colnames(temp_shared))){
    # Create a new column called sample_type if it is not already present
     temp_shared <- temp_shared %>% mutate(sampleID = sample) %>% 
        select(-sample) %>% select(sampleID, everything())
      
    } else{
      # reads in data (assumes tab delimited) if sample_type present
      temp_shared <- temp_shared %>% 
        filter(sample_type == "stool") %>% mutate(sampleID = sample) %>% 
        select(-sample) %>% select(sampleID, everything())
    }
    
  }
  # prints message to stdout updating on completion
  print(paste("Completed Reading in: ", i, ending, " data", sep = ""))
  # returns the necessary temp_shared file
  return(temp_shared)
}



# Function to grab groupings from DMM
grab_dmm_groups <- function(dataTable, metaData, kvalue = 2, seed_value = 1234567){
  # dataTable is the matrix to be analyzed
  # metaData is the respective meta data table for which grouping will be drawn from
  # kvalue represents the number of groups (components to make)
  # seed_value is the random number at which the seed will start
  
  # load and/or install needed package
  loadLibs("DirichletMultinomial")
  # run the test on the respective data
  dmm_test <- dmn(dataTable, k = kvalue, verbose = T, seed = seed_value)
  # pull the probablities for the first group
  dmm_groups <- dmm_test@group[, 1]
  # add column that assign groupings to the metadata table based on G1 probabilities
  metaData <- metaData %>% 
    mutate(dmm_groups = ifelse(dmm_groups >= 0.5, invisible("g1"), invisible("g2")))
  # return the modified metadata file with the group information added
  return(metaData)
  
}


# Function to create a 2x2 table and run a fisher test
get_fisher_pvalue <- function(metaData){
  
  options(show.error.messages = FALSE)
          
  summary_stat_value <- try(
    
   fisher.test(table(metaData$disease, metaData$dmm_groups))$p.value
  )
  
  options(show.error.messages = TRUE)
  
  summary_stat_value <- ifelse(is.numeric(summary_stat_value), 
                               invisible(summary_stat_value), invisible(NA))
  
  return(summary_stat_value)
  
}


# Control function to get all the data
get_data <- function(i){
  
  sample_names <- (get_file(i, "data/process/", "_genera_shared.csv") %>% 
                     mutate(sample_names = rownames(.)) %>% select(sample_names))[, "sample_names"]
  
  sub_genera_data <- get_file(i, "data/process/", "_subsample_genera.csv", rows_present = F, 
                              vec_of_rownames = sample_names) %>% 
    as.data.frame() %>% mutate(sample_ID = rownames(.)) %>% 
    select(sample_ID, everything())
  
  study_meta <- get_file(i, "data/process/", ".metadata", rows_present = F,  
                         "stool", metadata = T) %>% 
    mutate(disease = ifelse(disease == "polyp", invisible("control"), invisible(disease)))
  
  if(length(rownames(study_meta)) < length(rownames(sub_genera_data))){
    
    sub_genera_data <- sub_genera_data %>% slice(match(study_meta$sampleID, sample_ID))
    
  } else{
    
    study_meta <- study_meta %>% slice(match(sub_genera_data$sample_ID, sampleID))
  }
  
  sample_names <- sub_genera_data$sample_ID
  sub_genera_data <- as.matrix(sub_genera_data %>% select(-sample_ID))
  rownames(sub_genera_data) <- sample_names
  
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta)
  
  return(dataList)
  
}








##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


stool_study_data <- mapply(get_data, stool_sets, SIMPLIFY = F)

pvalues <- c()
cl <- makeCluster(2)
registerDoParallel(cl)



# Gets stool final sets
pvalues <- foreach(i=1:length(stool_sets)) %dopar% {
  
  library(dplyr)
  
  testData <- study_data[[stool_sets[i]]]
  
  study_meta <- grab_dmm_groups(apply(testData[["sub_genera_data"]], 2, function(x) round(x)), 
                                testData[["study_meta"]])
  
  pvalues <- rbind(pvalues, c(stool_sets[i], get_fisher_pvalue(study_meta)))
  
}

final_stats <- as.data.frame(pvalues[[5]], stringsAsFactors = F) %>% 
  bind_rows(as.data.frame(pvalues[[6]], stringsAsFactors = F)) %>% 
  rename(study = V1, pvalue = V2) %>% 
  mutate(pvalue = as.numeric(pvalue), bh = p.adjust(pvalue, method = "BH"))



# Write out the results of k = 2 dmm
write.csv(final_stats, "data/process/tables/stool_dmm_results.csv", 
          row.names = F)





















