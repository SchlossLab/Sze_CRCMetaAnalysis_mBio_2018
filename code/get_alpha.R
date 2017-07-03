### Generate Alpha Diversity Comparisons
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "car", "ggplot2"))

# Need to loop through and do this for every data set
datasets <- c("ahn", "baxter", "brim", "burns", "chen", "dejea", "flemer", "geng", 
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

# set up list to store the stool data
all_stool <- list(wang = c(), brim = c(), weir = c(), ahn = c(), zeller = c(), baxter = c())
all_st_metdata <- list(wang = c(), brim = c(), weir = c(), ahn = c(), zeller = c(), baxter = c())
all_st_z_alpha <- list(wang = c(), brim = c(), weir = c(), ahn = c(), zeller = c(), baxter = c())

# Run the for loop to process the data
for(i in 1:length(stool_sets)){

  # Load in alpha metrics to be used for stool
  all_stool[[i]] <- read.delim(paste("data/process/", stool_sets[i], "/", stool_sets[i], 
                                   ".groups.ave-std.summary", sep = ""), 
                             header = T, stringsAsFactors = F) %>% 
    filter(method == "ave") %>% 
    mutate(group = as.character(group)) %>% 
    select(group, sobs, shannon, shannoneven)
  
  # Load in metadata and match
  all_st_metdata[[i]] <- read.delim(
    paste("data/process/", 
          stool_sets[i], "/", stool_sets[i], ".metadata", sep = ""), 
    header = T, stringsAsFactors = F) %>% 
    mutate(sample = as.character(sample)) %>% 
    slice(match(all_stool[[i]]$group, sample))
  
  # Can use scale to z-score normalize
  all_st_z_alpha[[i]] <- all_stool[[i]] %>% 
    inner_join(all_st_metdata[[i]], by = c("group" = "sample")) %>% 
    mutate_at(vars(sobs:shannoneven), function(x) as.vector(scale(x))) %>% 
    mutate(set = stool_sets[i])
  
}

# Combine data together

stool_z_combined <- c()

for(i in 1:length(stool_sets)){
  
  stool_z_combined <- stool_z_combined %>% bind_rows(all_st_z_alpha[[i]])
}


# convert sex
stool_z_combined <- stool_z_combined %>% 
  mutate(sex = gsub("female", "f", sex), 
         sex = gsub("male", "m", sex), 
         sex = gsub("<not provided>", NA, sex))

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

