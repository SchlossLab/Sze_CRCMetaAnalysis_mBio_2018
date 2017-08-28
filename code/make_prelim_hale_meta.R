### Make Hale metadata pretty
### Re-organizing so that only the necessary columns are needed
### Marc Sze

# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("dplyr", "tidyr"))


# Load in metadata file
orig_meta <- read.delim("data/process/hale/hale_metadata.txt", 
                        stringsAsFactors = F, header = T)

# Modify table

new_meta <- orig_meta %>% 
  select(Description, Diagnosis1, gender, RACE, AGE, smoke.1.Y., PRIOR_CA, REL_CA, REL_PLYP) %>% 
  rename(sampleID = Description, sex = gender, age = AGE, 
         disease = Diagnosis1, race = RACE, smoker = smoke.1.Y.) %>% 
  mutate(white = ifelse(race == "White", invisible(1), 
                        ifelse(race == "Patient Refused" | race == "Unknown" | race == "blank", 
                               invisible(NA), invisible(0))), 
         disease = ifelse(grepl("ad", disease) == T, invisible("polyp"), invisible(disease)), 
         sex = ifelse(sex == "MALE", invisible("m"), 
                      ifelse(sex == 0, invisible("f"), invisible(NA))))

# Grab final variables to be used in analysis
good_meta <- new_meta %>% 
  select(sampleID, disease, sex, white, age)

# Write out the file
write.csv(good_meta, "data/process/hale/hale_metadata.csv", row.names = F)

