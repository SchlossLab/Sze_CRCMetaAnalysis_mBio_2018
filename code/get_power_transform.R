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
stool_sets <- c("wang", "brim", "weir", "ahn", "zeller", "baxter")

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
  
  if(!("sex" %in% colnames(all_metdata))){
    
    all_metdata <- all_metdata %>% mutate(sex = NA)
  }
  
  if(!("age" %in% colnames(all_metdata))){
    
    all_metdata <- all_metdata %>% mutate(age = NA)
  }
  
  if(!("bmi" %in% colnames(all_metdata))){
    
    all_metdata <- all_metdata %>% mutate(bmi = NA)
  }
  
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


# Function to calculate lambda
calc_lambda <- function(original_df, transformed_data){
  
  filtered_df <- original_df %>% 
    select(sobs, shannon, shannoneven)

  if(transformed_data[1, 1] > 0){
    
    lambda_values <- log(transformed_data[1, ])/log(filtered_df[1, ])
  } else{
    
    lambda_values <- log(-1*transformed_data[1, ])/log(filtered_df[1, ])
  }
  
  
  return(unname(t(lambda_values[1, ])[, 1]))
}



# Function to transform needed data, get the lambda and p-values for normality
get_lambda <- function(test_df){
  
  new_power_data <- test_df %>% 
    mutate_at(vars(sobs, shannon, shannoneven), 
              function(x) transformTukey(x, plotit = FALSE))
    
  test_new_data <- new_power_data %>% select(sobs, shannon, shannoneven)
  
  optimum_lambdas <- calc_lambda(test_df, test_new_data)
    
  pvalues <- apply(test_new_data, 2, function(x) shapiro.test(x)$p.value)
  
  values <- cbind(optimum_lambdas, pvalues)
  
  lambda_data <- list(summary_stats = values, transformed_data = new_power_data)
  
  return(lambda_data)
  
}


# Create combined stool data table
stool_data <- get_combined_table(c(stool_sets, both_sets), "stool")
transformed_stool <- lapply(stool_data, get_lambda)

# Create combined tissue data table
tissue_data <- get_combined_table(c(tissue_sets, both_sets), "tissue")
transformed_tissue <- lapply(tissue_data, get_lambda)

# solving for lambda
  # transformed_num = num ^ lambda
  # log(transformed_num) = lambda*log(num)
  # log(transformed_num) / log(num) = lambda








