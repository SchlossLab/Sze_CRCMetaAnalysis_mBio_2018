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
stool_sets <- c("wang", "weir", "ahn", "zeller", "baxter", "flemer", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
# Ignore chen for stool since there is only one case
both_sets <- c("chen", "flemer")

# The package DirichletMultinomial really messes with dplyr and makes it almost impossible to
# use things such as slice and match in combination. Make sure to do all table 
# transformations before loading and using it.



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
  # metaData should be the modified metadata file with the added groupings column
  
  # supress any errors for this analysis (some dmm groups will have NA)
  options(show.error.messages = FALSE)
  # tries to generate pvalue from fisher test otherwise outputs chr of error        
  summary_stat_value <- try(
    # compares proportions of case and cancer within the two groups  
   fisher.test(table(metaData$disease, metaData$dmm_groups))$p.value
  )
  # turns off suppression of error messages
  options(show.error.messages = TRUE)
  # adds the value or adds NA depending on whether a chr is present
  summary_stat_value <- ifelse(is.numeric(summary_stat_value), 
                               invisible(summary_stat_value), invisible(NA))
  # returns the pvalue
  return(summary_stat_value)
  
}


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
                         "stool", metadata = T) %>% 
    mutate(disease = ifelse(disease == "polyp", invisible("control"), invisible(disease)))
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
  sub_genera_data <- as.matrix(sub_genera_data %>% select(-sample_ID))
  rownames(sub_genera_data) <- sample_names
  # creates a list file with both data sets
  dataList <- list(sub_genera_data = sub_genera_data, 
                   study_meta = study_meta)
  # returns the combined list file
  return(dataList)
  
}


##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

# reads in all the stool data into one list
stool_study_data <- mapply(get_data, stool_sets, SIMPLIFY = F)

# assign needed values and processors for the analysis
pvalues <- c()
cl <- makeCluster(2)
registerDoParallel(cl)



# Gets stool final sets by running the analysis on 2 different processors
pvalues <- foreach(i=1:length(stool_sets)) %dopar% {
  
  library(dplyr)
  
  testData <- stool_study_data[[stool_sets[i]]]
  
  study_meta <- grab_dmm_groups(apply(testData[["sub_genera_data"]], 2, function(x) round(x)), 
                                testData[["study_meta"]])
  
  pvalues <- data_frame(study = stool_sets[i], pvalue = get_fisher_pvalue(study_meta))
  
}

# combines the seperate data together from the two processors and adds the bh correction
final_stats <- bind_rows(pvalues) %>% 
  mutate(bh = p.adjust(pvalue, method = "BH"))



# Write out the results of k = 2 dmm
write.csv(final_stats, "data/process/tables/stool_dmm_results.csv", 
          row.names = F)





















