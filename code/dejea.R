#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease

shared <- read.table("data/process/dejea/combined.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T)

metadata <- read.csv("data/process/dejea/dejea.seqdata.convscrconly.csv", stringsAsFactors=F)

#Get only relevant samples
rownames(shared) <- shared$Group
rownames(metadata) <- metadata$Run_s
shared <- shared[rownames(metadata), ]

#Check that everything is okay
stopifnot(shared$Group == metadata$Run_s)

sample <- metadata$Run_s
white <- metadata$white
disease <- metadata$disease
age <- metadata$age
sex <- metadata$sex

metadata <- cbind(sample=sample, white=white, disease=disease, age=age, sex=sex)

write.table(shared, file="data/process/dejea/dejea.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/dejea/dejea.metadata", quote=F, sep='\t', row.names=F)
 
