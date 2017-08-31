### Get the dmm groups - tissue
### Merge the group data with relevant metadata
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "vegan", "foreach", "doParallel"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
# Manually add geng where needed
tissue_sets <- c("dejea", "sana", "burns")

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



# Function to create a list of sampleIDs and disease
make_list <- function(i, dataTable_name){
  
  dataTable <- get(dataTable_name)
  
  new_table <- dataTable %>% filter(study == i) %>% 
    select(group, id, disease, sample_type)
  
  if(length(rownames(new_table)) != 0){
    
    return(new_table)
  }
  
}


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
  
  # returns the combined list file
  return(sub_genera_data)
  
}


# Function to match metadata with sample data
make_match <- function(i, dataList, metaList){
  
  dataTable <- dataList[[i]]
  metaTable <- metaList[[i]]
  
  if(length(rownames(metaTable)) < length(rownames(dataTable))){
    
    dataTable <- dataTable %>% slice(match(metaTable$group, sample_ID))
    
  } else {
    
    metaTable <- metaTable %>% slice(match(dataTable$sample_ID, group))
  }
  
  dataTable <- make_dmm_nice(dataTable)
  
  combined_data <- list(dataTable = dataTable, 
                        metaTable = metaTable)
  
  return(combined_data)
  
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


# Function to create rownames and make matrix
make_dmm_nice <- function(data_table){
  
  # re assigns the rown names while removing the extra column used for sorting
  sample_names <- data_table$sample_ID
  data_table <- as.matrix(data_table %>% select(-sample_ID))
  rownames(data_table) <- sample_names
  
  return(data_table)
}



##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################


# Read in specific data tables with samples that are matched and unmatched

tissue_matched <- read.csv("data/process/tables/alpha_tissue_matched_data.csv", 
                           header = T, stringsAsFactors = F) %>% 
  mutate(matchings = ifelse(disease == "cancer" | disease == "polyp", 1, 0))

tissue_unmatched <- read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                             header = T, stringsAsFactors = F) %>% 
  mutate(disease = gsub("adenoma", "polyp", disease))


matched_meta <- mapply(make_list, c("dejea", "burns", "geng"), "tissue_matched", SIMPLIFY = F)
unmatched_meta <- mapply(make_list, c(tissue_sets, both_sets), "tissue_unmatched", SIMPLIFY = F)


matched_data <- mapply(get_data, c("dejea", "burns", "geng"))
unmatched_data <- mapply(get_data, c(tissue_sets, both_sets))


matched_sets <- sapply(c("dejea", "burns", "geng"), 
                       function(x) make_match(x, matched_data, matched_meta), simplify = F)

unmatched_sets <- sapply(c(tissue_sets, both_sets), 
                       function(x) make_match(x, unmatched_data, unmatched_meta), simplify = F)


rm(tissue_matched, tissue_unmatched, matched_meta, matched_data, 
   unmatched_meta, unmatched_data)


# assign needed values and processors for the analysis
pvalues <- c()
matched_tissue <- c("dejea", "burns", "geng")
unmatched_tissue <- c(tissue_sets, "flemer")
cl <- makeCluster(2)
registerDoParallel(cl)



# Gets tissue matched final sets by running the analysis on 2 different processors
pvalues <- foreach(i=1:length(matched_tissue)) %dopar% {
  
  library(dplyr)
  
  testData <- matched_sets[[matched_tissue[i]]][["dataTable"]]
  metaData <- matched_sets[[matched_tissue[i]]][["metaTable"]]
  
  study_meta <- grab_dmm_groups(apply(testData, 2, function(x) round(x)), metaData)
  
  pvalues <- rbind(pvalues, c(matched_tissue[i], get_fisher_pvalue(study_meta)))
  
}

# combines the seperate data together from the two processors and adds the bh correction
matched_final_stats <- as.data.frame(pvalues[[2]], stringsAsFactors = F) %>% 
  bind_rows(as.data.frame(pvalues[[3]], stringsAsFactors = F)) %>% 
  rename(study = V1, pvalue = V2) %>% 
  mutate(pvalue = as.numeric(pvalue), bh = p.adjust(pvalue, method = "BH"))


rm(pvalues)
pvalues <- c()

# Gets tissue unmatched final sets by running the analysis on 2 different processors
pvalues <- foreach(i=1:length(unmatched_tissue)) %dopar% {
  
  library(dplyr)
  
  testData <- unmatched_sets[[unmatched_tissue[i]]][["dataTable"]]
  metaData <- unmatched_sets[[unmatched_tissue[i]]][["metaTable"]]
  
  study_meta <- grab_dmm_groups(apply(testData, 2, function(x) round(x)), metaData)
  
  pvalues <- rbind(pvalues, c(unmatched_tissue[i], get_fisher_pvalue(study_meta)))
  
}


# combines the seperate data together from the two processors and adds the bh correction
unmatched_final_stats <- as.data.frame(pvalues[[2]], stringsAsFactors = F) %>% 
  bind_rows(as.data.frame(pvalues[[4]], stringsAsFactors = F)) %>% 
  rename(study = V1, pvalue = V2) %>% 
  mutate(pvalue = as.numeric(pvalue), bh = p.adjust(pvalue, method = "BH"))



# Need to modify chen to eliminate rows without values

chen_meta <- unmatched_sets[["chen"]][["metaTable"]]
chen_data <- unmatched_sets[["chen"]][["dataTable"]]

rd_chen_data <- apply(chen_data, 2, function(x) round(x))



















