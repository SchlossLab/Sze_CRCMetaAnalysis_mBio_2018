# Generate oligo table to be used to create an oligo file

#Load needed libraries
library(dplyr)


test <- read.table("data/raw/wang/test.txt", stringsAsFactors = F)


total_cols <- length(test[1, ])

samples <- test[, seq(1, total_cols, by = 2)]
barcodes <- test[, seq(2, total_cols, by = 2)]

test <- cbind(t(samples), t(barcodes))

colnames(test) <- c("samples", "oligo")

test <- as.data.frame(test) %>% 
  mutate(barcodes = rep("BARCODE", length(sample))) %>% 
  select(barcodes, oligo, samples)

write.table(test, "data/process/wang/wang.oligos", sep = "\t", 
            quote = FALSE, col.names = F, row.names = F)
