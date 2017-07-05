#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease


# Load in needed functions and libraries
source('code/functions.R')

loadLibs(c("dplyr", "tidyr"))

shared <- read.table("data/process/dejea/combined.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared", 
                     header=T, stringsAsFactors = F)


metadata1 <- read.csv("data/process/dejea/dejea_metadata.csv", 
                      header = T, stringsAsFactors=F)

metadata2 <- read.delim("data/process/dejea/SraRunTable.txt", 
                        header = T, stringsAsFactors = F)

# Merge the metadata together and create needed columns
temp_data <- metadata2 %>% separate(Sample_Name_s, 
                               into = c("patient_id", "disease", "location")) %>% 
  select(Run_s, patient_id) %>% 
  inner_join(metadata1, by = "patient_id") %>% 
  mutate(white = ifelse(race == "caucasian" | race == "hispanic", 1, 0)) %>% 
  rename(sample = Run_s, disease = patient_type, site = tumor_site) %>% 
  select(-patient_id) %>% 
  mutate(matched = rep("y", length(sample)))


# Align shared file with metadata
temp_shared <- shared %>% slice(match(Group, temp_data$sample))


#Check that everything is okay
stopifnot(temp_shared$Group == temp_data$sample)

write.table(temp_shared, file="data/process/dejea/dejea.shared", quote=F, sep='\t', row.names=F)

write.table(temp_data, file="data/process/dejea/dejea.metadata", quote=F, sep='\t', row.names=F)
 
