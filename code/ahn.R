#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease
	# Ethinicity breakdown is as follows:
		# 1 = White not Hispanic
		# 2 = Black not Hispanic
		# 3 = Hispanic
		# 5 = Asian or Pacific Islander

shared <- read.table("data/process/ahn/combined.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", stringsAsFactors=F, header=T)

metadata <- read.csv("data/process/ahn/AhnData.csv", stringsAsFactors=F)

demodata <- read.table("data/process/ahn/phs000884.v1.pht004601.v1.p1.c1.Gut_Microbiome_Controls_Subject_Phenotypes.HMB-MDS.txt", skip=9, stringsAsFactors=F, header=T)

demodata <- demodata[match(metadata$gap_subject_id_s, demodata$dbGaP_Subject_ID), ]
demodata$white[which(demodata$RACE == 1)] <- 1
demodata$white[which(demodata$RACE > 1)] <- 0
shared <- shared[match(metadata$Run_s, shared$Group), ]

metadata$disease[grep("Yes", metadata$subject_is_affected_s)] <- "cancer"
metadata$disease[grep("No", metadata$subject_is_affected_s)] <- "control"

stopifnot(shared$Group == metadata$Run_s)

sample <- metadata$Run_s
white <- demodata$white
disease <- metadata$disease

metadata <- cbind(sample=sample, white=white, disease=disease)

write.table(shared, file="data/process/ahn/ahn.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/ahn/ahn.metadata", quote=F, sep='\t', row.names=F)
 
