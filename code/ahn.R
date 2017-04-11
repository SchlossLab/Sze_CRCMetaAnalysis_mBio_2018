#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease
	# Ethinicity breakdown is as follows:
		# 1 = White not Hispanic
		# 2 = Black not Hispanic
		# 3 = Hispanic
		# 5 = Asian or Pacific Islander

library("dplyr")
library("tidyr")

shared <- read.delim("data/process/ahn/combined.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared", 
                     stringsAsFactors=F, header=T)

metadata <- read.csv("data/process/ahn/AhnData.csv", stringsAsFactors=F)

demodata <- read.table("data/process/ahn/phs000884.v1.pht004601.v1.p1.c1.Gut_Microbiome_Controls_Subject_Phenotypes.HMB-MDS.txt", 
                       skip=10, stringsAsFactors=F, header=T) %>% 
  rename(gap_subject_id_s = dbGaP_Subject_ID)

#merge demo and meta into one table with variables of interest
combined_meta <- inner_join(demodata, metadata, by = "gap_subject_id_s") %>% 
  mutate(white = ifelse(RACE == 1, 1, 0), 
         disease = ifelse(subject_is_affected_s == "Yes", invisible("cancer"), invisible("control"))) %>% 
  select(Run_s, disease, AGE, sex_s, bmi) %>% 
  rename(sample = Run_s, age = AGE, sex = sex_s)


shared <- shared[match(combined_meta$sample, shared$Group), ]


stopifnot(shared$Group == combined_meta$sample)


write.table(shared, file="data/process/ahn/ahn.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/ahn/ahn.metadata", quote=F, sep='\t', row.names=F)
 
