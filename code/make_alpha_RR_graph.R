### Code to create needed alpha diversity RR graphs
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "wesanderson"))

# Load needed data tables (adenoma)
adn_all_stool <- read_csv("data/process/tables/alpha_adn_RR_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_adn_RR_ind_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V1-V3", 3), rep("V4", 6), rep("V3-V5", 3)))

adn_all_tissue <- read_csv("data/process/tables/alpha_adn_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_adn_RR_ind_tissue_results.csv"))

# Load in needed data tables (carcinoma)
crc_all_stool <- read_csv("data/process/tables/alpha_RR_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_RR_ind_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V3", 3), rep("V4", 3), rep("V3-V4", 3), 
                    rep("V4", 6), rep("V3-V5", 3), rep("V3-V4", 3)))

crc_all_tissue <- read_csv("data/process/tables/alpha_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_RR_ind_tissue_results.csv"))


#### zeller = V4 - #CDB38B
#### hale = V3-V5 - #FFA500
#### brim = V1-V3 - #8B5A00
#### baxter = V4 - #CDB38B
#### weir = V4 - #CDB38B
#### wang = V3 - #B8860B
### flemer = V3-V4 - #E3CF57
### ahn = V3-V4 - #E3CF57

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

adn_all_stool %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "zeller", "hale", "brim", "baxter"), 
                        labels = c( "Combined", "Zeller", "Hale", "Brim", "Baxter"))) %>% 
  filter(measure == "sobs") %>% 
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-2.2, 2.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + 
  scale_color_manual(values = c( '#000000', '#8B5A00', '#FFA500', '#CDB38B'))


crc_all_stool %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "zeller", "weir", "wang", "hale", 
                                   "flemer", "baxter", "ahn"), 
                        labels = c( "Combined", "Zeller", "Weir", "Wang",  "Hale", 
                                    "Flemer", "Baxter", "Ahn")), 
         region = factor(region, 
                         levels = c("combined", "V3", "V4", "V3-V4", "V3-V5"))) %>%  
  filter(measure == "sobs") %>% 
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-2.2, 2.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + 
  scale_color_manual(values = c( '#000000', '#B8860B', '#CDB38B', '#E3CF57', '#FFA500')) + 
  annotate("text", label = paste("Carcinoma"), x = -2.1, y = 8.4, size = 3)




##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################
