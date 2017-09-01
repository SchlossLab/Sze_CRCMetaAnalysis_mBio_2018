### Generate Alpha Diversity Comparisons
### With power transformation
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "rcompanion"))

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("wang", "brim", "weir", "ahn", "zeller", "baxter", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")


get_data <- function(i, sample_source){
  # i represents a character vector with data sets that should be worked through
  # sample_source represents the type of sample e.g. stool or tissue
  # splitting will be done based on this call 
  
  # Load in alpha metrics to be used for stool
  all_data <- read.delim(paste("data/process/", i, "/", i, 
                               ".groups.ave-std.summary", sep = ""), 
                         header = T, stringsAsFactors = F) %>% 
    filter(method == "ave") %>% 
    mutate(group = as.character(group)) %>% 
    select(group, sobs, shannon, shannoneven)
  
  
  # Load in metadata and match
  all_metdata <- read.delim(
    paste("data/process/", i, "/", i, ".metadata", sep = ""), 
    header = T, stringsAsFactors = F) 
  
  # Create a new column called sample_type if it is not already present
  if(!("sample_type" %in% colnames(all_metdata))){
    
    all_metdata <- all_metdata %>% mutate(sample_type = sample_source)
  }
  # Create a new column called sex if it is not already present
  if(!("sex" %in% colnames(all_metdata))){
    
    all_metdata <- all_metdata %>% mutate(sex = NA)
  }
  # Create a new column called age if it is not already present
  if(!("age" %in% colnames(all_metdata))){
    
    all_metdata <- all_metdata %>% mutate(age = NA)
  }
  # Create a new column called bmi if it is not already present
  if(!("bmi" %in% colnames(all_metdata))){
    
    all_metdata <- all_metdata %>% mutate(bmi = NA)
  }
  # Create a new column called white if it is not already present
  if(!("white" %in% colnames(all_metdata))){
    
    all_metdata <- all_metdata %>% mutate(white = NA)
  }
  
  # filter based on sample type and match metadata with all alpha data
  combined_data <- all_metdata %>%   
    mutate(sample = as.character(sample)) %>% 
    filter(sample_type == sample_source) %>% 
    slice(match(all_data$group, sample)) %>% 
    mutate(study = i) %>% 
    inner_join(all_data, by = c("sample" = "group"))
  
  # convert sex column so that all entries are uniform
  combined_df <- combined_data %>% 
    mutate(sex = gsub("f.*", "f", sex, ignore.case = T), 
           sex = gsub("m.*", "m", sex, ignore.case = T), 
           sex = gsub("^(?!m|f).*$", NA, sex, perl = T, ignore.case = T))
  
  # convert disease column so that all normal are control
  # and all adenoma are polyp
  combined_df <- combined_data %>% 
    mutate(disease = ifelse(disease == "normal", invisible("control"), 
                            ifelse(disease == "adenoma", invisible("polyp"), 
                                   invisible(disease))))
  
    # Select specific columns and rows for the final data table
  combined_df <- combined_df %>% 
    select(sample, sobs, shannon, shannoneven, disease, 
           white, sample_type, sex, age, bmi, study) %>% 
    filter(!is.na(disease)) %>% 
    rename(group = sample)
  

  return(combined_df)
}


# Create function to wrangle the data
get_combined_table <- function(datasets, sample_source){
  # i represents a character vector with data sets that should be worked through
  # sample_source represents the type of sample e.g. stool or tissue
  # splitting will be done based on this call 
  
  # Get alpha metrics with relevant metadata
  combined_data <- mapply(get_data, datasets, sample_source, SIMPLIFY = FALSE)
  
  # combine all the different data sets together
  #combined_df <- bind_rows(lapply(combined_data, function(x) as.data.frame(x)))
  
  # Write out completed data table
  return(combined_data)
  
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


# write out tables 
write_out_tables <- function(dataset, data_list, sample_type, table_type,  
                             directory = "data/process/tables/", further_sep = NULL){
  
  if(is.null(further_sep)){
    
    write.csv(data_list, 
              paste(directory, dataset, "_", sample_type, "_", table_type, ".csv", sep = ""), 
              row.names = F)
    
  } else{
    
    write.csv(data_list[[further_sep]], 
              paste(directory, dataset, "_", sample_type, "_", 
                    further_sep, "_", table_type, ".csv", sep = ""), 
              row.names = F)
  }

}



# Create combined stool data table
stool_data <- get_combined_table(c(stool_sets, both_sets), "stool")
transformed_stool <- lapply(stool_data, get_lambda)

# Write out stool_data
mapply(write_out_tables, c(stool_sets, both_sets), stool_data, "stool", "alpha_raw_values")
mapply(write_out_tables, c(stool_sets, both_sets), transformed_stool, 
       "stool", "alpha_raw_values", further_sep = "transformed_data")

write.csv(bind_rows(lapply(transformed_stool, function(x) 
  x[["summary_stats"]] %>% as.data.frame() %>% 
    mutate(alpha_measure = c("sobs", "shannon", "shannoneven"), 
           study = unique(x[["transformed_data"]]$study)) %>% 
    rename(shapiro_wilk_pvalue = pvalues))), 
  "data/process/tables/stool_power_transformation_summary.csv", row.names = F)



mapply(write_out_tables, c(stool_sets, both_sets), transformed_stool, 
       "stool", "alpha_raw_values", further_sep = "summary_stats")

# Create combined tissue data table
tissue_data <- get_combined_table(c(tissue_sets, both_sets), "tissue")
transformed_tissue <- lapply(tissue_data, get_lambda)

# Write out tissue_data
mapply(write_out_tables, c(tissue_sets, both_sets), tissue_data, "tissue", "alpha_raw_values")
mapply(write_out_tables, c(tissue_sets, both_sets), transformed_tissue, 
       "tissue", "alpha_raw_values", further_sep = "transformed_data")
mapply(write_out_tables, c(tissue_sets, both_sets), transformed_tissue, 
       "tissue", "alpha_raw_values", further_sep = "summary_stats")

# solving for lambda
  # transformed_num = num ^ lambda
  # log(transformed_num) = lambda*log(num)
  # log(transformed_num) / log(num) = lambda








