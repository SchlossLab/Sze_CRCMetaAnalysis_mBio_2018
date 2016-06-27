#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease

shared <- read.table("data/process/lu/lu.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared", header=T)

metadata <- read.csv("data/process/lu/luData.csv", stringsAsFactors=F)

rownames(shared) <- shared[, "Group"]
rownames(metadata) <- metadata[, "Run_s"]
matched_metadata <- metadata[rownames(shared), ]

stopifnot(shared$Group == matched_metadata$Run_s)

unmatchedsamples <- matched_metadata[which(
  matched_metadata$divider == "C" | matched_metadata$divider == "A"), ] 
umatchedshared <- shared[rownames(unmatchedsamples), ]

stopifnot(umatchedshared$Group == unmatchedsamples$Run_s)

matchedsamples <- matched_metadata[which(
  matched_metadata$divider == "B" | matched_metadata$divider == "A"), ]
matchedshared <- shared[rownames(matchedsamples), ]

stopifnot(matchedshared$Group == matchedsamples$Run_s)

sample <- matchedsamples$Run_s
white <- rep(0,length(rownames(matchedsamples)))
disease <- matchedsamples$disease
pairs <- matchedsamples$Sample_Name_s

matchedmetadata <- cbind(sample=sample, white=white, disease=disease, pairs=pairs)

write.table(matchedshared, file="data/process/lu/lu.matched.shared", quote=F, sep='\t', row.names=F)

write.table(matchedmetadata, file="data/process/lu/lu.matched.metadata", quote=F, sep='\t', row.names=F)


sample <- unmatchedsamples$Run_s
white <- rep(0,length(rownames(unmatchedsamples)))
disease <- unmatchedsamples$disease
pairs <- unmatchedsamples$Sample_Name_s

unmatchedmetadata <- cbind(sample=sample, white=white, disease=disease, pairs=pairs)

write.table(umatchedshared, file="data/process/lu/lu.unmatched.shared", quote=F, sep='\t', row.names=F)

write.table(unmatchedmetadata, file="data/process/lu/lu.unmatched.metadata", quote=F, sep='\t', row.names=F)

 
