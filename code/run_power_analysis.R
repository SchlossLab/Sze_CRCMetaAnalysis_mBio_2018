### Power Analysis for Stool and tissue data
### Get power for crc and adn for the different study sets with respect to effect size
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
tissue_sets <- c("dejea", "lu", "geng", "sana", "burns")



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
  # joins data together based on needed sampleID column and then removes it keeping on disease
  sub_genera_data <- study_meta %>% 
    inner_join(sub_genera_data, by = c("sampleID" = "sample_ID")) %>% 
    select(-sampleID)
  
  # returns the combined list file
  return(sub_genera_data$disease)
  
}

# Function to generate needed proportions dependent on disease severity
generate_ES_values <- function(i, dataList, diseaseType){
  # i represents the study
  # dataList is the list of interest from the read in of data counts
  # diseaseType is the disease of interest (e.g. "cancer" or "polyp")
  
  # generate a temporary vector based on i
  tempData <- dataList[[i]]
  # Check to see if i is any of these studies
  if(i %in% c("brim", "lu", "flemer_t") & diseaseType == "cancer"){
    # assign null to the read out since they do not exist for that disease type
    values <- NULL
  } else{
    # Check what disease type we are looking for
    if(diseaseType == "cancer"){
      # change polyps to control for cancer comparison
      tempData <- gsub("polyp", "control", tempData)
      # generate counts for each group
      overall_counts <- as.data.frame(table(tempData))
      # assign the needed info from the table to be read out
      values <- c(
        p.dis = filter(overall_counts, tempData == "cancer")[, "Freq"] / sum(overall_counts$Freq), 
        n_case = filter(overall_counts, tempData == "cancer")[, "Freq"], 
        n_control = filter(overall_counts, tempData == "control")[, "Freq"])
      
    } else{
      # Check if i is any of these data sets 
      if(i %in% c("brim", "zeller", "baxter", "hale", "lu", "flemer_t")){
        # Remove all cancer samples
        tempData <- tempData[tempData != "cancer"]
        # generate the counts for each group
        overall_counts <- as.data.frame(table(tempData))
        # assign the needed info from the table to be read out
        values <- c(
          p.dis = filter(overall_counts, tempData == "polyp")[, "Freq"] / sum(overall_counts$Freq), 
          n_case = filter(overall_counts, tempData == "polyp")[, "Freq"], 
          n_control = filter(overall_counts, tempData == "control")[, "Freq"])
        
      } else{
        # if it lands here it has neither cancer or polyp present in the data
        values <- NULL
      }
    }
  }
  # Return the values needed for pwr analysis to the global environment
  return(values)
}


# Function to get the actual power and predicted n
generate_pwr <- function(i, dataList){
  # i is the study of interest
  # dataList is the list of interest containing the proportions and counts
  
  # grab the vector of data of interest from list
  tempData <- dataList[[i]]
  # The % difference we are interested in
  differences <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
  # check to see if the study is null
  if(is.null(tempData)){
    # If it is null assign null to the values
    values <- NULL
  } else{
    # Run the proportion test for each % difference of interest and get the associated power
    study_power <- sapply(differences, 
                          function(x) 
                            pwr.2p2n.test(
                              h=ES.h(tempData["p.dis"], tempData["p.dis"] + x),
                                 n1=tempData["n_case"], n2=tempData["n_control"],
                                 alternative="two.sided")$power)
    # assign the % differences and column names
    names(study_power) <- differences
    
    # Generate the hypothetical balanced n needed for each % difference 
    needed_n <- sapply(differences, 
                       function(x) 
                         pwr.2p.test(
                           h=ES.h(tempData["p.dis"], tempData["p.dis"] + x), power=0.80)$n)
    # assign the percent difference as column names
    names(needed_n) <- differences
    # Create a final nice data table with all the necessary results
    values <- tibble(
      study_power = study_power, 
      pwr80_needed_n = needed_n, 
      effect_size = differences, 
      study = i)
  }
  # read out the results to the global environment
  return(values)
  
}




##############################################################################################
############### Run the actual programs to get the data (ALL Data) ###########################
##############################################################################################

# reads in all the stool and tissue data
stool_study_data <- mapply(get_data, stool_sets, "stool", SIMPLIFY = F)
tissue_study_data <- mapply(get_data, tissue_sets, "tissue", SIMPLIFY = F)

# The meta data for the flemer study tissue component is not in the typical metadata.csv file
flemer_tissue <- list(flemer_t = (read.csv("data/process/tables/alpha_tissue_unmatched_data.csv", 
                          header = T, stringsAsFactors = F) %>% 
  filter(study == "flemer"))[, "disease"])

# Combine all the data sets together since for this it doesn't matter if it is tissue or stool
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










