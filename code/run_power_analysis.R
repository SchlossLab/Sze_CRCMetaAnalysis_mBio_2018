### Power Analysis for Stool
### Get power for crc and adn for the different stool sets with respect to effect size
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')


# Load needed libraries
loadLibs(c("tidyverse", "pwr", "statmod"))

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
# Ignore brim since it only has polyps
stool_sets <- c("brim", "wang", "weir", "ahn", "zeller", "baxter", "hale", "flemer")

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
# Remove Lu since it only has polyps and no cancer cases
# Have to remove dejea since there are no controls
# Have to omit Chen since their is only 1 case
tissue_sets <- c("lu", "geng", "sana", "burns")



##############################################################################################
############### List of functions needed to run the analysis #################################
##############################################################################################

# Control function to get all the data, basically runs the above functions in a
# contained location withouth having to repeat them
get_data <- function(i, datatype){
  # i is the study of interest
  
  # gets original sample names
  sample_names <- (get_file(i, "data/process/", "_genera_shared.csv") %>% 
                     mutate(sample_names = rownames(.)) %>% select(sample_names))[, "sample_names"]
  # grabs subsampled data and assigns rownames from sample names to table
  sub_genera_data <- get_file(i, "data/process/", "_subsample_genera.csv", rows_present = F, 
                              vec_of_rownames = sample_names) %>% 
    as.data.frame() %>% mutate(sample_ID = rownames(.)) %>% 
    select(sample_ID, everything()) %>% mutate(sample_ID = as.character(sample_ID))
  # grabs the meta data and transforms polyp to control (polyp/control vs cancer) 
  study_meta <- get_file(i, "data/process/", ".metadata", rows_present = F,  
                         datatype, metadata = T) %>% 
    mutate(sampleID = as.character(sampleID)) %>% 
    select(sampleID, disease) %>% 
    mutate(disease = stringr::str_replace(disease, "normal", "control"), 
           disease = stringr::str_replace(disease, "adenoma", "polyp"))
  
  sub_genera_data <- study_meta %>% 
    inner_join(sub_genera_data, by = c("sampleID" = "sample_ID")) %>% 
    select(-sampleID)
  
  # returns the combined list file
  return(sub_genera_data$disease)
  
}

# Function to generate needed proportions dependent on disease severity
generate_ES_values <- function(i, dataList, diseaseType){
  
  tempData <- dataList[[i]]
  
  if(i %in% c("brim", "lu", "flemer_t") & diseaseType == "cancer"){
    
    values <- NULL
  } else{
    
    if(diseaseType == "cancer"){
      
      tempData <- gsub("polyp", "control", tempData)
      overall_counts <- as.data.frame(table(tempData))
      values <- c(
        p.dis = filter(overall_counts, tempData == "cancer")[, "Freq"] / sum(overall_counts$Freq), 
        n_case = filter(overall_counts, tempData == "cancer")[, "Freq"], 
        n_control = filter(overall_counts, tempData == "control")[, "Freq"])
      
    } else{
      
      if(i %in% c("brim", "zeller", "baxter", "hale", "lu", "flemer_t")){
        
        tempData <- tempData[tempData != "cancer"]
        overall_counts <- as.data.frame(table(tempData))
        values <- c(
          p.dis = filter(overall_counts, tempData == "polyp")[, "Freq"] / sum(overall_counts$Freq), 
          n_case = filter(overall_counts, tempData == "polyp")[, "Freq"], 
          n_control = filter(overall_counts, tempData == "control")[, "Freq"])
        
      } else{
        
        values <- NULL
      }
    }
  }
  
  return(values)
}


# Function to get the actual power and predicted n
generate_pwr <- function(i, dataList){
  
  tempData <- dataList[[i]]
  differences <- c(0.01, 0.05, 0.10, 0.15, 0.20)
  
  if(is.null(tempData)){
    
    values <- NULL
  } else{
    
    study_power <- sapply(differences, 
                          function(x) 
                            pwr.2p2n.test(
                              h=ES.h(tempData["p.dis"], tempData["p.dis"] + x),
                                 n1=tempData["n_case"], n2=tempData["n_control"],
                                 alternative="two.sided")$power)
    
    names(study_power) <- differences
    
    # Assumes a balanced design
    needed_n <- sapply(differences, 
                       function(x) 
                         pwr.2p.test(
                           h=ES.h(tempData["p.dis"], tempData["p.dis"] + x), power=0.80)$n)
    
    names(needed_n) <- differences
   
    values <- tibble(
      study_power = study_power, 
      pwr80_needed_n = needed_n, 
      effect_size = differences, 
      study = i)
  }
  
  return(values)
  
}




##############################################################################################
############### Run the actual programs to get the data (ALL Data) ###########################
##############################################################################################

# reads in all the stool data into one list
stool_study_data <- mapply(get_data, stool_sets, "stool", SIMPLIFY = F)

tissue_study_data <- mapply(get_data, tissue_sets, "tissue", SIMPLIFY = F)

flemer_tissue <- list(flemer_t = (read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                          header = T, stringsAsFactors = F) %>% 
  filter(study == "flemer"))[, "disease"])


combined_data <- c(stool_study_data, tissue_study_data, flemer_tissue)

# Generate needed info for effect size testing
cancer_info <- sapply(c(stool_sets, tissue_sets, "flemer_t"), 
                      function(x) generate_ES_values(x, combined_data, diseaseType = "cancer"))

polyp_info <- sapply(c(stool_sets, tissue_sets, "flemer_t"), 
                     function(x) generate_ES_values(x, combined_data, diseaseType = "polyp"))

# Generate the needed power measurements 
cancer_pwr_data <- sapply(c(stool_sets, tissue_sets, "flemer_t"), 
                          function(x) generate_pwr(x, cancer_info)) %>% 
  bind_rows()

polyp_pwr_data <- sapply(c(stool_sets, tissue_sets, "flemer_t"), 
                          function(x) generate_pwr(x, polyp_info)) %>% 
  bind_rows()

# Write out the needed data
write.csv(cancer_pwr_data, file="data/process/tables/cancer_predicted_pwr_and_n.csv", 
          row.names=F)
write.csv(polyp_pwr_data, file="data/process/tables/adn_predicted_pwr_and_n.csv", 
          row.names=F)










