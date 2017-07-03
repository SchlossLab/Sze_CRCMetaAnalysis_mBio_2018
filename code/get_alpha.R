### Generate Alpha Diversity Comparisons
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "car"))

# Need to loop through and do this for every data set


# Load in alpha metrics to be used
test_data <- read.delim("data/process/ahn/ahn.groups.ave-std.summary", 
                      header = T, stringsAsFactors = F) %>% 
  filter(method == "ave") %>% 
  select(group, sobs, shannon, shannoneven)

# Load in metadata and match
meta_data <- read.delim("data/process/ahn/ahn.metadata", 
                        header = T, stringsAsFactors = F) %>% 
  slice(match(test_data$group, sample))

# Can use scale to z-score normalize
test_data <- test_data %>% inner_join(meta_data, by = c("group" = "sample")) %>% 
  mutate_at(vars(sobs:shannoneven), function(x) as.vector(scale(x)))
  





# Can use powerTransform instead to try and correct out skew

