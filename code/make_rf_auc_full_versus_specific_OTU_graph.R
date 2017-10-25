### Code to measure AUCs of full data set and select OTU dataset
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Load needed data tables (adenoma)
adn_tissue <- read_csv("data/process/tables/adn_tissue_rf_otu_random_comparison_summary.csv") %>% 
  mutate(type = c("matched", "unmatched"), model_type = "full") %>% 
  bind_rows(read_csv("data/process/tables/adn_tissue_rf_select_otu_random_comparison_summary.csv") %>% 
              mutate(type = c("matched", "unmatched"), model_type = "select")) %>% 
  gather(key = measure_type, value = AUC, act_mean_auc, act_sd_auc, rand_mean_auc, rand_sd_auc)

adn_all_stool <- read_csv("data/process/tables/adn_stool_rf_otu_random_comparison_summary.csv") %>% 
  mutate(model_type = "full") %>% 
  bind_rows(read_csv("data/process/tables/adn_stool_rf_select_otu_random_comparison_summary.csv") %>% 
              mutate(model_type = "select")) %>% 
  gather(key = measure_type, value = AUC, act_mean_auc, act_sd_auc, rand_mean_auc, rand_sd_auc)

# Load in needed data tables (carcinoma)
crc_tissue_matched <- read_csv("data/process/tables/matched_tissue_rf_otu_random_comparison_summary.csv") %>% 
  mutate(type = "matched", model_type = "full") %>% 
  bind_rows(read_csv("data/process/tables/matched_tissue_rf_select_otu_random_comparison_summary.csv") %>% 
              mutate(type = "matched", model_type = "select")) %>% 
  gather(key = measure_type, value = AUC, act_mean_auc, act_sd_auc, rand_mean_auc, rand_sd_auc)

crc_tissue_unmatched <- 
  read_csv("data/process/tables/unmatched_tissue_rf_otu_random_comparison_summary.csv") %>% 
  mutate(type = "unmatched", model_type = "full") %>% 
  bind_rows(read_csv("data/process/tables/unmatched_tissue_rf_select_otu_random_comparison_summary.csv") %>% 
              mutate(type = "unmatched", model_type = "select")) %>% 
  gather(key = measure_type, value = AUC, act_mean_auc, act_sd_auc, rand_mean_auc, rand_sd_auc)

crc_all_stool <- read_csv("data/process/tables/stool_rf_otu_random_comparison_summary.csv") %>% 
  mutate(model_type = "full") %>% 
  bind_rows(read_csv("data/process/tables/stool_rf_select_otu_random_comparison_summary.csv") %>% 
              mutate(model_type = "select")) %>% 
  gather(key = measure_type, value = AUC, act_mean_auc, act_sd_auc, rand_mean_auc, rand_sd_auc)



##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

### Study Colors by Viridis 
##  library(scales) 
## show_col(viridis_pal()(16))
# flemer - #440154FF
# lu - #FDE725FF
# burns - #453581FF
# chen - #3D4D8AFF
# sana - #1F998AFF
# dejea - #2B748EFF
# geng - #CBE11EFF
# brim - #34618DFF
# zeller - #FDE725FF
# baxter - #481D6FFF
# hale - #67CC5CFF
# wang - #97D83FFF
# weir - #24878EFF
# ahn - #40BC72FF

adn_tissue %>% 
  mutate(type = factor(type, 
                       levels = c("unmatched", "matched"), 
                       labels = c("Unmatched Tissue", "Matched Tissue")), 
         model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("flemer", "lu"), 
                        labels = c("Flemer", "Lu"))) %>% 
  filter(grepl("sd", measure_type) != T, grepl("rand", measure_type) != T) %>% 
  ggplot(aes(model_type, AUC, color = study, group = study)) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) +
  geom_point(size = 3.5, show.legend = F) + 
  geom_line(show.legend = F) + facet_grid(. ~ type) + coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#440154FF', '#FDE725FF')) + 
  annotate("text", label = paste("Adenoma"), x = 0.6, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))



##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

