#Output shared file that has the same rows as the metadata file. 
#Requirements:
	#  *Rows must be in the same order
	#  *Metadata must contain sample id, white, and disease

# Necessary Extra information
  # ERR260325 or TCA1 is for healthy control stool
  # ERR260326 or TCA2 is for CRC stool
  # this information was used to craft the raw.metad.csv file

# Email Primary contact:
#  Tiffany Weir at Tiffany.Weir@colostate.edu

shared <- read.table("data/process/weir/combined.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared", 
                     stringsAsFactors=F, header=T)

metadata <- read.csv("data/process/weir/weir.raw.metadata.csv", 
                     stringsAsFactors=F, header = T)

stopifnot(shared$Group == metadata$SampleID)

sample <- metadata$SampleID
white <- metadata$white
disease <- metadata$disease
sex <- metadata$sex
age <- metadata$age

metadata <- cbind(sample=sample, white=white, disease=disease, sex=sex, age=age)

write.table(shared, file="data/process/weir/weir.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/weir/weir.metadata", quote=F, sep='\t', row.names=F)
 
