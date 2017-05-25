# Generate oligo table to be used to create an oligo file

library("dplyr")
library("tidyr")

test <- read.delim("data/raw/chen/run1.txt", header = F, stringsAsFactors = F) %>% 
  separate(V2, c("sample", "barcode")) %>% 
  select(sample, barcode) %>% filter(is.na(barcode) != TRUE) %>% 
  rename(oligo = barcode) %>% mutate(barcode = rep("BARCODE", length(rownames(.)))) %>% 
  select(barcode, oligo, sample)

test2 <- read.delim("data/raw/chen/run2.txt", header = F, stringsAsFactors = F) %>% 
  separate(V2, c("sample", "barcode")) %>% 
  select(sample, barcode) %>% filter(is.na(barcode) != TRUE) %>% 
  rename(oligo = barcode) %>% mutate(barcode = rep("BARCODE", length(rownames(.)))) %>% 
  select(barcode, oligo, sample)


write.table(test, "data/process/chen/run1.oligos", sep = "\t", 
            quote = FALSE, col.names = F, row.names = F)
write.table(test2, "data/process/chen/run2.oligos", sep = "\t", 
          quote = FALSE, col.names = F, row.names = F)

