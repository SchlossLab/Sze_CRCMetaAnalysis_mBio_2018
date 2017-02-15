#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease

shared <- read.table(
	"data/process/zeller/zeller.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared", 
	header=T, stringsAsFactors=F)

metadata <- read.csv("data/raw/Zeller/ZellerMetadata.csv", header=T, stringsAsFactors=F)

#Get only relevant samples
rownames(shared) <- shared$Group
rownames(metadata) <- metadata$Run_s
metadata <- metadata[rownames(shared), ]

metadata$dx[which(metadata$Diagnosis_s == "Normal")] <- "control" 
metadata$dx[which(metadata$Diagnosis_s == "Small adenoma" | 
	metadata$Diagnosis_s == "Large adenoma")] <- "polyp"
metadata$dx[which(metadata$Diagnosis_s == "Cancer")] <- "cancer"

#Check that everything is okay
stopifnot(shared$Group == metadata$Run_s)

sample <- metadata$Run_s
disease <- metadata$dx
age <- metadata$age_s
sex <- metadata$sex_s

metadata <- cbind(sample=sample, disease=disease, age=age, sex=sex)

write.table(shared, file="data/process/zeller/zeller.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/zeller/zeller.metadata", quote=F, sep='\t', row.names=F)
 
