#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("tidyverse"))

shared <- read_tsv("data/process/dejea/combined.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared")


metadata1 <- read_csv("data/process/dejea/dejea_metadata.csv")

metadata2 <- read_tsv("data/process/dejea/SraRunTable.txt")

# Merge the metadata together and create needed columns
temp_data <- metadata2 %>% separate(Sample_Name_s, 
                               into = c("patient_id", "disease", "location")) %>% 
  select(Run_s, patient_id, disease, location) %>% 
  inner_join(metadata1, by = "patient_id") %>% 
  mutate(white = ifelse(race == "caucasian" | race == "hispanic", 1, 0)) %>% 
  rename(sample = Run_s, site = tumor_site) %>% 
  select(-patient_type) %>% 
  
  mutate(disease = stringr::str_replace(disease, "Normal", "control"), 
         disease = stringr::str_replace(disease, "Tumor", "cancer"), 
         disease = stringr::str_replace(disease, "Polyp", "polyp"), 
    matched = rep("y", length(sample)))


# Align shared file with metadata
temp_shared <- shared %>% slice(match(Group, temp_data$sample))


#Check that everything is okay
stopifnot(temp_shared$Group == temp_data$sample)

write.table(temp_shared, file="data/process/dejea/dejea.shared", quote=F, sep='\t', row.names=F)

write_tsv(temp_data, "data/process/dejea/dejea.metadata")
 
