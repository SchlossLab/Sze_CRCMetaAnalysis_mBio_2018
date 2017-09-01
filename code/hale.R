### create metadata file for Hale, et al study 2016
### Marc Sze

### Some metadata is pulled from the study


# Load needed libraries
library(dplyr)
library(tidyr)

# Load needed data
metadata <- read.csv("data/process/hale/hale_metadata.csv", 
                     header = T, stringsAsFactors = F)

# This is too big to load on laptop
shared <- read.delim("data/process/hale/hale.trim.contigs.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared", 
                     header = T, stringsAsFactors = F)

# Match shared with tempData
select_meta <- metadata %>% slice(match(shared$Group, sampleID))


# Check order is okay
stopifnot(shared$Group == select_meta$samples)

select_meta <- select_meta %>% rename(sample = sampleID)

# Write out to table
write.table(shared, file="data/process/hale/hale.shared", quote=F, sep='\t', row.names=F)

write.table(select_meta, file="data/process/hale/hale.metadata", quote=F, sep='\t', row.names=F)

