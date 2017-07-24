### Generate Alpha Diversity Comparisons
### Z transform power transformed data
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr", "car", "ggplot2", "lme4"))

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


# Function to read in the respective transformed data tables for each study

get_transformed_data <- function(i, sampleType){
  
  data_list <- read.csv(paste("data/process/tables/", i, "_", sampleType, "_",  
                                   "transformed_data_alpha_raw_values.csv", sep = ""), 
                             header = T, stringsAsFactors = F)

  return(data_list)
}


# Function to apply a z-score normalization and and a new column with these values
# for sobs, shannon, and shannoneven

zscore_transform <- function(dataList, 
                             variables = c("sobs", "shannon", "shannoneven")){
  
  z_trans_data <- dataList %>% 
    mutate_at(vars(variables), funs(as.numeric(scale(.))))
  
  return(z_trans_data)
}


# Read in the respective data
stool_transformed_data <- mapply(get_transformed_data, 
                                 c(stool_sets, both_sets), "stool", SIMPLIFY = F)

# Apply z-score normalization to every data set
stool_ztrans_pwrtrans_data <- lapply(stool_transformed_data, zscore_transform)

# create a combined data table
combined_data <- bind_rows(lapply(stool_ztrans_pwrtrans_data, 
                                  function(x) as.data.frame(x) %>% 
                           mutate(group = as.character(group))))

# Test the stool data

t.test(filter(combined_data, disease == "cancer")[, "sobs"], 
       filter(combined_data, disease != "cancer")[, "sobs"], 
       alternative = "two.sided")

t.test(filter(combined_data, disease == "cancer")[, "shannon"], 
       filter(combined_data, disease != "cancer")[, "shannon"], 
       alternative = "two.sided")

t.test(filter(combined_data, disease == "cancer")[, "shannoneven"], 
       filter(combined_data, disease != "cancer")[, "shannoneven"], 
       alternative = "two.sided")

# Test if there is a difference based on adenoma being seperated from controls
TukeyHSD(aov(lm(combined_data$sobs ~ factor(combined_data$disease))))


# Account for data sets in analysis using a linear mixed-effect model
null_model <- lmer(combined_data$shannon ~ (1|combined_data$study), REML=FALSE)
disease_model <- lmer(combined_data$shannon ~ (1|combined_data$disease) + (1|combined_data$study), 
                      REML=FALSE)
anova(null_model, disease_model)[["Pr(>Chisq)"]][2]




