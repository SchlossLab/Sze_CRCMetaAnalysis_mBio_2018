### Generate Alpha Diversity Comparisons
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2"))

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


get_Z_data <- function(i, sample_source){
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
  
  # filter based on sample type and match metadata with all alpha data
  all_metdata <- all_metdata %>%   
    mutate(sample = as.character(sample)) %>% 
    filter(sample_type == sample_source) %>% 
    slice(match(all_data$group, sample))
  
  # Can use scale to z-score normalize
  all_z_alpha <- all_data %>% 
    inner_join(all_metdata, by = c("group" = "sample")) %>% 
    mutate_at(vars(sobs, shannon, shannoneven), function(x) as.vector(scale(x))) %>% 
    mutate(study = i)
  
  return(all_z_alpha)
}


# Create function to wrangle the data
get_combined_table <- function(datasets, sample_source){
  # i represents a character vector with data sets that should be worked through
  # sample_source represents the type of sample e.g. stool or tissue
    # splitting will be done based on this call 
  
  # Get Z-score transformed alpha with relevant metadata
  all_z_alpha <- mapply(get_Z_data, datasets, sample_source)
  
  # combine all the different data sets together
  data_z_combined <- bind_rows(lapply(all_z_alpha, function(x) as.data.frame(x)))

  
  # convert sex column so that all entries are uniform
  data_z_combined <- data_z_combined %>% 
    mutate(sex = gsub("f.*", "f", sex, ignore.case = T), 
           sex = gsub("m.*", "m", sex, ignore.case = T), 
           sex = gsub("^(?!m|f).*$", NA, sex, perl = T, ignore.case = T))
  
  # Write out completed data table
  return(data_z_combined)
  
}

# Create function to get histograms of each measure
generate_histogram_alpha <- function(data_set, sample_type){
  
  # Set up the variables that will be used
  alpha_measures <- c("sobs", "shannon", "shannoneven")
  alpha_names <- c("Richness", "Shannon Diversity", "Evenness")
  overall_hist <- list(sobs = c(), shannon = c(), shannoneven = c())
  
  # Visualize the distributions for the stool data set
  for(i in 1:length(alpha_measures)){
    
    # Generates the actual plot
    overall_hist[[i]] <- ggplot(data_set, aes(data_set[, alpha_measures[i]])) + 
      geom_histogram(color = "black", fill = "white", bins = 40) + theme_bw() + 
      ylab("Counts") + xlab(paste("Z-score Normalized ", alpha_names[i], sep = "")) + 
      scale_y_continuous(limits = c(0, 110), expand = c(0, 0))
    
    # Save graphs in exploratory notebook
    ggsave(paste("exploratory/notebook/", 
                 sample_type, "_", alpha_measures[i], "_combined_histogram.tiff", sep = ""), 
           overall_hist[[i]], width=6, height = 7, dpi = 300)
  }
  
  return(overall_hist)
}


# Create combined stool data table
stool_data <- get_combined_table(stool_sets, "stool") %>% 
  bind_rows(get_combined_table(both_sets, "stool")) %>% 
  select(group, sobs, shannon, shannoneven, disease, white, sample_type, sex, 
         age, bmi, study) %>% 
  filter(!is.na(disease))


# Create combined tissue data table
tissue_data <- get_combined_table(tissue_sets, "tissue") %>% 
  bind_rows(get_combined_table(both_sets, "tissue")) %>% 
  select(group, sobs, shannon, shannoneven, white, disease, matched, sample_type, 
         age, sex, site, stage, size_mm, bmi, study)


# Create graphs
stool_graphs <- generate_histogram_alpha(stool_data, sample_type = "stool")
tissue_graphs <- generate_histogram_alpha(tissue_data, sample_type = "tissue")


# Write out combined data
write.csv(stool_data, "data/process/tables/stool_zscore_alpha_combined.csv", row.names = F)
write.csv(tissue_data, "data/process/tables/tissue_zscore_alpha_combined.csv", row.names = F)

# Note to self: 
# Can use powerTransform instead to try and correct out skew

