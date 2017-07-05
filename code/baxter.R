#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease

shared <- read.table(
	"data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.shared", 
	header=T)

metadata <- read.delim("data/process/baxter/metadata.tsv", header=T, stringsAsFactors=F)

#Get only relevant samples
rownames(shared) <- shared$Group
rownames(metadata) <- metadata$sample
shared <- shared[rownames(metadata), ]

metadata$dx[which(metadata$dx == "normal")] <- "control" 
metadata$dx[which(metadata$dx == "adenoma")] <- "polyp"

#Check that everything is okay
stopifnot(shared$Group == metadata$sample)

sample <- metadata$sample
white <- metadata$White
disease <- metadata$dx
age <- metadata$Age
sex <- metadata$Gender
bmi <- metadata$BMI

metadata <- cbind(sample=sample, white=white, disease=disease, age=age, sex=sex, bmi=bmi)

write.table(shared, file="data/process/baxter/baxter.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/baxter/baxter.metadata", quote=F, sep='\t', row.names=F)
 
