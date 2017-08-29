### Get the dmm groups 
### Merge the group data with relevant metadata
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "vegan"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("wang", "brim", "weir", "ahn", "zeller", "baxter")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")



# Function to read in needed genera file
get_file <- function(i, path_to_file, ending, rows_present=T, name_of_rows=1, 
                     vec_of_rownames = NULL, make_matrix = F, 
                     sample_source, metadata = F){
  
  if(metadata == F){
    
    if(rows_present == T){
      
      temp_shared <- read.csv(paste(path_to_file, i, "/", i, ending, sep = ""), 
                              header = T, stringsAsFactors = F, row.names = name_of_rows)
    } else{
      
      temp_shared <- read.csv(paste(path_to_file, i, "/", i, ending, sep = ""), 
                              header = T, stringsAsFactors = F)
      rownames(temp_shared) <- vec_of_rownames
    }
    
  } else{
    
    temp_shared <- read.delim(paste(path_to_file, i, "/", i, ending, sep = ""), 
                              header = T, stringsAsFactors = F)
    
    # Create a new column called sample_type if it is not already present
    if(!("sample_type" %in% colnames(temp_shared))){
      
     temp_shared <- temp_shared %>% mutate(sampleID = sample) %>% 
        select(-sample) %>% select(sampleID, everything())
      
    } else{
      
      temp_shared <- temp_shared %>% 
        filter(sample_type == "stool") %>% mutate(sampleID = sample) %>% 
        select(-sample) %>% select(sampleID, everything())
    }
    
  }
  
  print(paste("Completed Reading in: ", i, ending, " data", sep = ""))
  
  return(temp_shared)
}



# Function to grab groupings from DMM
grab_dmm_groups <- function(dataTable, metaData, kvalue = 2, seed_value = 1234567){
  
  dmm_test <- dmn(dataTable, k = kvalue, verbose = T, seed = seed_value)
  
  dmm_groups <- dmm_test@group[, 1]
  
  metaData <- metaData %>% 
    mutate(dmm_groups = ifelse(dmm_groups >= 0.5, invisible("g1"), invisible("g2")))
  
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




### TO DO LIST ###

# filter data sets to make sure only stool or tissue being compared
# modify disease category so it is only control versus cancer
# remove polyp only data sets
# create a control function to read in data
# use lapply instead of the for loop





##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

pvalues <- c()

for(i in c("weir")){
  
  sample_names <- (get_file(i, "data/process/", "_genera_shared.csv") %>% 
    mutate(sample_names = rownames(.)) %>% select(sample_names))[, "sample_names"]
  
  sub_genera_data <- get_file(i, "data/process/", "_subsample_genera.csv", rows_present = F, 
             vec_of_rownames = sample_names) %>% 
    as.data.frame() %>% mutate(sample_ID = rownames(.)) %>% 
    select(sample_ID, everything())
  
  study_meta <- get_file(i, "data/process/", ".metadata", rows_present = F,  
                         "stool", metadata = T)
  
  if(length(rownames(study_meta)) < length(rownames(sub_genera_data))){
    
    sub_genera_data <- sub_genera_data %>% slice(match(study_meta$sampleID, sample_ID))
    
  } else{
    
    study_meta <- study_meta %>% slice(match(sub_genera_data$sample_ID, sampleID))
  }

  
  #study_meta <- grab_dmm_groups(sub_genera_data, study_meta)
  
  #pvalues <- c(pvalues, get_fisher_pvalue(study_meta))
  
}


#final_stats <- cbind(pvalue = pvalues, bh = p.adjust(pvalues, method = "BH"))
  





















