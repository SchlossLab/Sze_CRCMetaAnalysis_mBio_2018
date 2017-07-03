### create metadata file for wang, et al study 2012
### Marc Sze

### Some data is pulled from the study
### Other components pulled from https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP005150

# Load needed libraries
library(dplyr)

# Load needed data
samples <- read.delim("data/process/wang/wang.oligos", header = F)

shared <- read.delim("data/process/wang/combined.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared", 
                     header = T, stringsAsFactors = F)

# Check order is okay
stopifnot(shared$Group == samples$V3)

# Add metadata

metadata <- samples %>% 
  rename(sample = V3) %>% 
  mutate(disease = ifelse(grepl("C", sample) == TRUE, invisible("cancer"), invisible("control")), 
         white = rep(0, length(rownames(.)))) %>% 
  select(sample, disease, white)

# Write out to table

write.table(shared, file="data/process/wang/wang.shared", quote=F, sep='\t', row.names=F)

write.table(metadata, file="data/process/wang/wang.metadata", quote=F, sep='\t', row.names=F)


