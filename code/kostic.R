#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease
	# Ethinicity breakdown is as follows:
		# 1 = White not Hispanic
		# 2 = Black not Hispanic
		# 3 = Hispanic
		# 5 = Asian or Pacific Islander

library("tidyverse")

shared <- read_tsv("data/process/kostic/combined.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared")

metadata <- read_tsv("data/process/kostic/SraRunTable.txt")

#merge demo and meta into one table with variables of interest
good_meta <- metadata %>%  
  mutate(white = NA,  
         disease = ifelse(isolation_source_s == "Tumor", invisible("cancer"), invisible("control")), 
         bmi = NA, 
         age = NA) %>% 
  select(Run_s, Sample_name_s, white, disease, age, host_sex_s, bmi) %>% 
  rename(sample = Run_s, sex = host_sex_s)


good_meta <- good_meta %>% slice(match(shared$Group, sample))


stopifnot(shared$Group == good_meta$sample)


write_tsv(shared, "data/process/kostic/kostic.shared")

write_tsv(good_meta, "data/process/kostic/kostic.metadata")
 
