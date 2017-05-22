#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease

shared <- read.table("data/process/lu/lu.trim.contigs.unique.good.filter.unique.pick.pick.opti_mcc.unique_list.shared", 
                     stringsAsFactors = F, header=T)

metadata <- read.csv("data/process/lu/luData.csv", stringsAsFactors=F, header = T)

rownames(shared) <- shared[, "Group"]
rownames(metadata) <- metadata[, "Run_s"]
matched_metadata <- metadata[rownames(shared), ]

stopifnot(shared$Group == matched_metadata$Run_s)

sample <- matched_metadata$Run_s
white <- rep(0,length(rownames(matched_metadata)))
disease <- matched_metadata$disease
matched <- matched_metadata$matched

matchedmetadata <- cbind(sample=sample, white=white, disease=disease, matched=matched)

write.table(shared, file="data/process/lu/lu.shared", quote=F, sep='\t', row.names=F)

write.table(matchedmetadata, file="data/process/lu/lu.metadata", quote=F, sep='\t', row.names=F)



