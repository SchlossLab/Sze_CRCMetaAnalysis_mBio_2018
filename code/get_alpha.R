### Generate Alpha Diversity Comparisons
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "car"))


test_data <- read.delim("data/process/ahn/ahn.groups.ave-std.summary", 
                      header = T, stringsAsFactors = F) %>% 
  filter(method == "ave")



# Can use scale to z-score normalize


# Can use powerTransform instead to try and correct out skew

