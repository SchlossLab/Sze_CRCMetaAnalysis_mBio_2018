### Create Modified Metagenome Download List
### Marc Sze

library(dplyr)

data_set <- read.csv("data/raw/zeller_metagenome_metadata.csv", header = T, stringsAsFactors = F)

test <- data_set %>% filter(InsertSize_l == 250) %>% distinct(host_subject_id_s, .keep_all = TRUE)

write.csv(test, "data/raw/metagenome_samples_to_keep_metadata.csv", row.names = F)
