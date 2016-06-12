#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease

shared <- read.table("data/brim/test/combined.trim.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T)

metadata <- read.csv("data/brim/brim.metadata.csv", stringAsFactors=F)

stopifnot(shared$Group == metadata$sample)

sample <- metadata$sample
white <- metadata$white
disease <- metadata$disease

metadata <- cbind(sample=sample, white=white, disease=disease)

write.table(shared, file="data/brim/test/brim.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/brim/test/brim.metadata", quote=F, sep='\t', row.names=F)
 
