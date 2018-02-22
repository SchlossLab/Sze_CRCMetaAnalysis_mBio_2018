### Code to create alpha graphs of control vs carcinoma
### Stool specific measures only
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Stool Only sets
# Hale, Wang, Brim, Weir, Ahn, Zeller, Baxter
stool_sets <- c("wang", "brim", "weir", "ahn", "zeller", "baxter", "hale")

# Both Tissue and Stool
# flemer sampletype = biopsy or stool
# chen sample_type = tissue or stool
# Flemer, Chen
both_sets <- c("flemer", "chen")





##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

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

# Create combined stool data table
stool_data <- get_combined_table(c(stool_sets, both_sets), "stool") %>% 
  bind_rows() %>% 
  select(-sobs, -white, -sample_type, -sex, -age, -bmi) %>% 
  filter(disease != "polyp") %>% 
  filter(study != "brim", study != "chen") %>% 
  mutate(study = factor(study, 
                        levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                        labels = c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller")), 
         disease = factor(disease, 
                          levels = c("control", "cancer"), 
                          labels = c("Control", "Carcinoma")))

adn_stool_data <- get_combined_table(c(stool_sets, both_sets), "stool") %>% 
  bind_rows() %>% 
  select(-sobs, -white, -sample_type, -sex, -age, -bmi) %>% 
  filter(study == "brim" | study == "zeller" | study == "baxter" | study == "hale") %>% 
  filter(disease != "cancer") %>% 
  mutate(study = factor(study, 
                        levels = c("baxter", "brim", "hale", "zeller"), 
                        labels = c("Baxter", "Brim", "Hale", "Zeller")), 
         disease = factor(disease, 
                          levels = c("control", "polyp"), 
                          labels = c("Control", "Adenoma")))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


adn_even <- adn_stool_data %>% 
  ggplot(aes(study, shannoneven, color = disease, group = disease)) + 
  geom_point(position = position_dodge(width = 0.6), alpha = 0.5, show.legend = T) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", position = position_dodge(width = 0.6), 
               size = 0.5, width = 0.4) +
  labs(x = "", y = "Evenness") + 
  scale_color_manual(name = "", 
                     values = c('#4169E1', '#DC143C')) + 
  theme_bw() + coord_flip(ylim = c(0, 1)) + ggtitle("A") + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8), 
        legend.background = element_rect(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))




even <- stool_data %>% 
  ggplot(aes(study, shannoneven, color = disease, group = disease)) + 
  geom_point(position = position_dodge(width = 0.6), alpha = 0.5, show.legend = T) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", position = position_dodge(width = 0.6), 
               size = 0.5, width = 0.4) +
  labs(x = "", y = "Evenness") + 
  scale_color_manual(name = "", 
                     values = c('#4169E1', '#DC143C')) + 
  theme_bw() + coord_flip(ylim = c(0, 1)) + ggtitle("B") + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8), 
        legend.background = element_rect(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))




shannon_d <- stool_data %>% 
  ggplot(aes(study, shannon, color = disease, group = disease)) + 
  geom_point(position = position_dodge(width = 0.6), alpha = 0.5) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", position = position_dodge(width = 0.6), 
               size = 0.5, width = 0.4) +
  labs(x = "", y = "Shannon Diversity") + 
  scale_color_manual(name = "", 
                     values = c('#4169E1', '#DC143C')) + 
  theme_bw() + coord_flip(ylim = c(0, 7)) + ggtitle("C") + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8), 
        legend.background = element_rect(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

combined_graph <- grid.arrange(adn_even, even, shannon_d, ncol = 3, nrow = 1)

ggsave("results/figures/Figure1.pdf", combined_graph, width = 10.5, height = 4, dpi = 300)








