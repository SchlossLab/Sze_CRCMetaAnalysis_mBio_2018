### create metadata file for Flemer, et al study 2016
### Marc Sze

### Some metadata is pulled from the study


# Load needed libraries
library(dplyr)
library(tidyr)

# Load needed data
metadata <- read.csv("data/process/flemer/gutjnl-2015-309595-inline-supplementary-material-1.csv", 
                     header = T, stringsAsFactors = F, skip = 2)

seq_match_data <- read.csv("data/process/flemer/flemer_2016_gut.csv", 
                           header = T, stringsAsFactors = F)

# Use this as a proxy for groups and filtering.
#proxy_shared <- read.delim2("flemer.trim.contigs.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.count.summary", 
#                           header = F, stringsAsFactors = F) %>% 
#  rename(Group = V1, seq_total = V2) 


# This is too big to load on laptop
shared <- read.delim("data/process/flemer/flemer.trim.contigs.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared", 
                     header = T, stringsAsFactors = F) %>% mutate(Group = seq_match_data$ID)


# Convert periods to _
shared <- shared %>% mutate(Group = gsub("\\.", "_", Group))
metadata <- metadata %>% mutate(ID = gsub("\\.", "_", ID))

# Match shared with tempData
select_meta <- metadata %>% slice(match(shared$Group, ID))

# Check order is okay
stopifnot(shared$Group == select_meta$samples)

# Write out to table

write.table(shared, file="data/process/flemer/flemer.shared", quote=F, sep='\t', row.names=F)

write.table(select_meta, file="data/process/flemer/flemer.metadata", quote=F, sep='\t', row.names=F)


