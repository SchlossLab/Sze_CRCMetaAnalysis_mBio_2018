### create metadata file for wang, et al study 2012
### Marc Sze

### Some data is pulled from the study
### Other components pulled from https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP005150

# Load needed libraries
library(dplyr)
library(tidyr)

# Some warings are thrown but these are not important
# Turn the warnings off
options(warn = -1)

# Load needed data
test <- read.delim("data/raw/chen/run1.txt", header = F, stringsAsFactors = F) %>% 
  separate(V2, c("sample", "barcode"))

test2 <- read.delim("data/raw/chen/run2.txt", header = F, stringsAsFactors = F) %>% 
  separate(V2, c("sample", "barcode")) 


shared <- read.delim("data/process/chen/chen.shared", 
                     header = T, stringsAsFactors = F)


# need to create seperate disease vectors
groups1 <- test$V1[2:length(test$V1)]
sample_id1 <- test$sample[1:(length(test$sample)-1)] # removes one blank value at the end

groups2 <- test2$V1[2:length(test2$V1)]
sample_id2 <- test2$sample[1:(length(test2$sample)-1)] # removes one blank value at the end

# recombine to create data tables
temp1 <- as.data.frame(cbind(sample_id1, groups1)) %>% 
  rename(sample = sample_id1, classes = groups1) %>% 
  mutate(sample_type = ifelse(grepl("stool", classes) == TRUE, invisible("stool"), invisible("tissue")), 
         disease = ifelse(grepl("colorectal cancer", classes) == TRUE, invisible("cancer"), invisible("control")), 
         white = rep(0, length(rownames(.)))) %>% 
  select(sample, sample_type, disease, white)

temp2 <- as.data.frame(cbind(sample_id2, groups2)) %>% 
  rename(sample = sample_id2, classes = groups2) %>% 
  mutate(sample_type = ifelse(grepl("swab", classes) == TRUE, invisible("swab"), invisible("tissue")), 
         disease = ifelse(grepl("colorectal cancer", classes) == TRUE, invisible("cancer"), invisible("control")), 
         white = rep(0, length(rownames(.)))) %>% 
  select(sample, sample_type, disease, white)


tempData <- temp1 %>% bind_rows(temp2)

# Match shared with tempData
select_meta <- tempData %>% slice(match(shared$Group, sample)) %>% 
  mutate(sex = rep(NA, length(sample)))

# Check order is okay
stopifnot(shared$Group == select_meta$sample)

# Write out to table

write.table(shared, file="data/process/chen/chen.shared", quote=F, sep='\t', row.names=F)

write.table(select_meta, file="data/process/chen/chen.metadata", quote=F, sep='\t', row.names=F)


