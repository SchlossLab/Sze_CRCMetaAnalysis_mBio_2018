### Generate Alpha Diversity Comparisons
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "car", "ggplot2"))

# Need to loop through and do this for every data set
total_datasets <- c("ahn", "baxter", "brim", "burns", "chen", "dejea", "flemer", "geng", 
              "lu", "sana", "wang", "weir", "zeller")

# Tissue Only sets
# Lu, Dejea, Sana, Burns, Geng
tissue_sets <- c("lu", "dejea", "sana", "burns", "geng")

# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("wang", "brim", "weir", "ahn", "zeller", "baxter")

# Both Tissue and Stool
# Flemer, Chen
both_sets <- c("flemer", "chen")

# Create function to wrangle the data
get_combined_table <- function(datasets, sample_source){
  
  # set up list to store the stool data
  all_data <- as.list(datasets)
  names(all_data) <- datasets
  all_metdata <- all_data
  all_z_alpha <- all_data
  
  # Run the for loop to process the data
  for(i in 1:length(datasets)){
    
    # Load in alpha metrics to be used for stool
    all_data[[i]] <- read.delim(paste("data/process/", datasets[i], "/", datasets[i], 
                                       ".groups.ave-std.summary", sep = ""), 
                                 header = T, stringsAsFactors = F) %>% 
      filter(method == "ave") %>% 
      mutate(group = as.character(group)) %>% 
      select(group, sobs, shannon, shannoneven)
    
    # Load in metadata and match
    all_metdata[[i]] <- read.delim(
      paste("data/process/", 
            datasets[i], "/", datasets[i], ".metadata", sep = ""), 
      header = T, stringsAsFactors = F) %>% 
      mutate(sample = as.character(sample)) %>% 
      slice(match(all_data[[i]]$group, sample))
    
    # Can use scale to z-score normalize
    all_z_alpha[[i]] <- all_data[[i]] %>% 
      inner_join(all_metdata[[i]], by = c("group" = "sample")) %>% 
      mutate_at(vars(sobs:shannoneven), function(x) as.vector(scale(x))) %>% 
      mutate(set = datasets[i], sample_type = sample_source)
    
  }
  
  # Combine data together
  data_z_combined <- c()
  
  for(i in 1:length(datasets)){
    
    print(datasets[i])
    data_z_combined <- data_z_combined %>% 
      bind_rows(mutate(all_z_alpha[[i]]))
  }
  
  
  # convert sex
  data_z_combined <- data_z_combined %>% 
    mutate(sex = gsub("female", "f", sex), 
           sex = gsub("male", "m", sex), 
           sex = gsub("<not provided>", NA, sex))
  
  return(data_z_combined)
  
}


stool_data <- get_combined_table(stool_sets, "stool")
tissue_data <- get_combined_table(tissue_sets, "tissue")







# Set up the variables that will be used
alpha_measures <- c("sobs", "shannon", "shannoneven")
alpha_names <- c("Richness", "Shannon Diversity", "Evenness")
overall_hist <- list(sobs = c(), shannon = c(), shannoneven = c())

# Visualize the distributions
for(i in 1:length(alpha_measures)){
  
  overall_hist[[i]] <- ggplot(stool_z_combined, aes(stool_z_combined[, alpha_measures[i]])) + 
    geom_histogram(color = "black", fill = "white", bins = 40) + theme_bw() + 
    ylab("Counts") + xlab(paste("Z-score Normalized ", alpha_names[i], sep = "")) + 
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0))
  
  # Save graphs in exploratory notebook
  ggsave(paste("exploratory/notebook/", alpha_measures[i], "_combined_histogram.tiff", sep = ""), 
         overall_hist[[i]], width=6, height = 7, dpi = 300)
}


# Write out combined data
write.csv(stool_z_combined, "data/process/tables/stool_zscore_alpha_combined.csv", row.names = F)


# Can use powerTransform instead to try and correct out skew

