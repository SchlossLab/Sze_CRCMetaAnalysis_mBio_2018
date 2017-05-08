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

shared <- read.delim("data/process/sanapareddy/combined.unique.good.filter.unique.precluster.pick.precluster.pick.pick.opti_mcc.unique_list.shared", 
                     stringsAsFactors=F, header=T)

metadata <- read.delim("data/process/sanapareddy/sequence_group_data.txt", stringsAsFactors=F, header = T) %>% 
  mutate(disease = ifelse(caseContol == "control", invisible("control"), invisible("cancer"))) %>% 
  select(sample, disease)


shared <- shared[match(metadata$sample, shared$Group), ]


stopifnot(shared$Group == metadata$sample)


write.table(shared, file="data/process/sanapareddy/sana.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/sanapareddy/sana.metadata", quote=F, sep='\t', row.names=F)
 
