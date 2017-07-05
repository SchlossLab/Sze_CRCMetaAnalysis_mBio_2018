#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease

# Load files
shared <- read.table("data/process/burns/burns.trim.contigs.good.unique.good.filter.unique.pick.pick.opti_mcc.unique_list.shared", 
                     header=T, stringsAsFactors = F)

metadata1 <- read.csv("data/process/burns/burnsMetadata.csv", stringsAsFactors=F)

metadata2 <- read.csv("data/process/burns/burnsMetadata2.csv", header = T, stringsAsFactors = F, 
                      nrows = 88)

# Load Libraries
library(dplyr)
library(tidyr)


#Split and merge columns two metadata files
metadata <- rename(metadata1, SampleID = Sample_Name_s) %>% 
  inner_join(metadata2, by = "SampleID") %>% 
  mutate(disease = ifelse(Description == "tumor", invisible("cancer"), invisible("control"))) %>% 
  slice(match(shared$Group, Run_s)) %>% 
  select(Run_s, disease, age_s, sex_s, Site, Stage) %>% 
  rename(sample = Run_s, age = age_s, sex = sex_s, site = Site, stage = Stage) %>% 
  mutate(matched = rep("y", length(sample)), stage = as.numeric(stage))
 
# Test to make sure everythin is matched
stopifnot(shared$Group == metadata$Run_s)

write.table(shared, file="data/process/burns/burns.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/burns/burns.metadata", quote=F, sep='\t', row.names=F)
 
