#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease

shared <- read.table("data/process/geng/combined.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T)

metadata <- read.csv("data/process/geng/GengData.csv", stringsAsFactors=F)

metadata$disease[grep("H", metadata$Sample_Name_s)] <- "control"
metadata$disease[grep("C", metadata$Sample_Name_s)] <- "cancer"

stopifnot(shared$Group == metadata$Run_s)

sample <- metadata$Run_s
white <- rep(0,length(rownames(metadata)))
disease <- metadata$disease

metadata <- cbind(sample=sample, white=white, disease=disease)

write.table(shared, file="data/process/geng/geng.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/geng/geng.metadata", quote=F, sep='\t', row.names=F)
 
