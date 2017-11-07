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
  gather(key = measure_type, value = AUC, act_mean_auc, rand_mean_auc, act_sd_auc, rand_sd_auc) %>% 
  separate(measure_type, c("model", "measure_type", "X1")) %>% 
  select(-X1) %>% 
  spread(key = measure_type, value = AUC) %>% 
  rename(AUC = mean, deviation = sd)

adn_all_stool <- read_csv("data/process/tables/adn_stool_rf_otu_random_comparison_summary.csv") %>% 
  mutate(model_type = "full") %>% 
  bind_rows(read_csv("data/process/tables/adn_stool_rf_select_otu_random_comparison_summary.csv") %>% 
              mutate(model_type = "select")) %>% 
  gather(key = measure_type, value = AUC, act_mean_auc, rand_mean_auc, act_sd_auc, rand_sd_auc) %>% 
  separate(measure_type, c("model", "measure_type", "X1")) %>% 
  select(-X1) %>% 
  spread(key = measure_type, value = AUC) %>% 
  rename(AUC = mean, deviation = sd)

# Load in needed data tables (carcinoma)
crc_tissue_matched <- read_csv("data/process/tables/matched_tissue_rf_otu_random_comparison_summary.csv") %>% 
  mutate(type = "matched", model_type = "full") %>% 
  bind_rows(read_csv("data/process/tables/matched_tissue_rf_select_otu_random_comparison_summary.csv") %>% 
              mutate(type = "matched", model_type = "select")) %>% 
  gather(key = measure_type, value = AUC, act_mean_auc, rand_mean_auc, act_sd_auc, rand_sd_auc) %>% 
  separate(measure_type, c("model", "measure_type", "X1")) %>% 
  select(-X1) %>% 
  spread(key = measure_type, value = AUC) %>% 
  rename(AUC = mean, deviation = sd) %>% 
  unite(combined_study, study, type, sep = "_")

crc_tissue_unmatched <- 
  read_csv("data/process/tables/unmatched_tissue_rf_otu_random_comparison_summary.csv") %>% 
  mutate(type = "unmatched", model_type = "full") %>% 
  bind_rows(read_csv("data/process/tables/unmatched_tissue_rf_select_otu_random_comparison_summary.csv") %>% 
              mutate(type = "unmatched", model_type = "select")) %>% 
  gather(key = measure_type, value = AUC, act_mean_auc, rand_mean_auc, act_sd_auc, rand_sd_auc) %>% 
  separate(measure_type, c("model", "measure_type", "X1")) %>% 
  select(-X1) %>% 
  spread(key = measure_type, value = AUC) %>% 
  rename(AUC = mean, deviation = sd) %>% 
  unite(combined_study, study, type, sep = "_")

crc_all_tissue <- crc_tissue_matched %>% bind_rows(crc_tissue_unmatched)


crc_all_stool <- read_csv("data/process/tables/stool_rf_otu_random_comparison_summary.csv") %>% 
  mutate(model_type = "full") %>% 
  bind_rows(read_csv("data/process/tables/stool_rf_select_otu_random_comparison_summary.csv") %>% 
              mutate(model_type = "select")) %>% 
  gather(key = measure_type, value = AUC, act_mean_auc, rand_mean_auc, act_sd_auc, rand_sd_auc) %>% 
  separate(measure_type, c("model", "measure_type", "X1")) %>% 
  select(-X1) %>% 
  spread(key = measure_type, value = AUC) %>% 
  rename(AUC = mean, deviation = sd)



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

adn_tissue_graph <- adn_tissue %>% 
  mutate(type = factor(type, 
                       levels = c("unmatched", "matched"), 
                       labels = c("Unmatched Tissue", "Matched Tissue")), 
         model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "CRC Associated\nGenera Community Only")), 
         study = factor(study, 
                        levels = c("flemer", "lu"), 
                        labels = c("Flemer", "Lu\n(matched)"))) %>% 
  filter(grepl("rand", model) != T) %>% 
  ggplot(aes(study, AUC, color = study, group = study)) + 
  geom_point(size = 3.5, show.legend = F) + 
  geom_errorbar(aes(ymin = AUC-deviation, ymax = AUC+deviation), 
                width = 0.25, size = 0.7, show.legend = F) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  facet_grid(. ~ model_type) + coord_flip(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#440154FF', '#FDE725FF')) + 
  annotate("text", label = paste("Adenoma"), x = 2.5, y = 0.05, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


adn_stool_graph <- adn_all_stool %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "CRC Associated\nGenera Community Only")), 
         study = factor(study, 
                        levels = c("baxter", "brim", "hale", "zeller"), 
                        labels = c("Baxter", "Brim", "Hale", "Zeller"))) %>% 
  filter(grepl("rand", model) != T) %>% 
  ggplot(aes(study, AUC, color = study, group = study)) + 
  geom_point(size = 3.5, show.legend = F) + 
  geom_errorbar(aes(ymin = AUC-deviation, ymax = AUC+deviation), 
                width = 0.25, size = 0.7, show.legend = F) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_flip(ylim = c(0, 1.05)) + 
  facet_grid(. ~ model_type) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#481D6FFF', '#1F998AFF', '#67CC5CFF', '#FDE725FF')) + 
  annotate("text", label = paste("Adenoma (Stool)"), x = 4.5, y = 0.1, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))



crc_tissue_graph <- crc_all_tissue %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "CRC Associated\nGenera Community Only")), 
         combined_study = factor(combined_study, 
                        levels = c("sana_unmatched", "geng_matched", "flemer_unmatched", 
                                   "dejea_matched", "chen_unmatched", "burns_unmatched",
                                   "burns_matched"), 
                        labels = c("Sanapareddy", "Geng\n(matched)", "Flemer", "Dejea\n(matched)", 
                                   "Chen", "Burns\n(matched)", "Burns"))) %>% 
  filter(grepl("rand", model) != T) %>% 
  ggplot(aes(combined_study, AUC, color = combined_study, group = combined_study)) + 
  geom_point(size = 3.5, show.legend = F) + 
  geom_errorbar(aes(ymin = AUC-deviation, ymax = AUC+deviation), 
                width = 0.25, size = 0.7, show.legend = F) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  facet_grid(. ~ model_type) + coord_flip(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#1F998AFF', '#CBE11EFF', '#440154FF', '#2B748EFF', 
                                '#3D4D8AFF', '#453581FF', '#453581FF')) + 
  annotate("text", label = paste("Carcinoma"), x = 7.35, y = 0.06, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.15, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))



crc_stool_graph <- crc_all_stool %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "CRC Associated\nGenera Community Only")), 
         study = factor(study, 
                        levels = c("zeller", "weir", "wang", "hale", "flemer", "baxter", "ahn"), 
                        labels = c("Zeller", "Weir", "Wang", "Hale", "Flemer", "Baxter", "Ahn"))) %>% 
  filter(grepl("rand", model) != T) %>% 
  ggplot(aes(study, AUC, color = study, group = study)) + 
  geom_point(size = 3.5, show.legend = F) + 
  geom_errorbar(aes(ymin = AUC-deviation, ymax = AUC+deviation), 
                width = 0.25, size = 0.7, show.legend = F) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_flip(ylim = c(0, 1.05)) + 
  facet_grid(. ~ model_type) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#FDE725FF', '#24878EFF', '#97D83FFF', '#67CC5CFF', 
                                '#440154FF', '#481D6FFF', '#40BC72FF')) + 
  annotate("text", label = paste("Carcinoma (Stool)"), x = 7.4, y = 0.1, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))



##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


tissue_auc_graph <- grid.arrange(adn_tissue_graph, crc_tissue_graph)


stool_auc_graph <- grid.arrange(adn_stool_graph, crc_stool_graph)

ggsave("results/figures/Figure5.pdf", 
       tissue_auc_graph, width = 7, height = 6, dpi = 300)

ggsave("results/figures/Figure4.pdf", 
       stool_auc_graph, width = 7, height = 6, dpi = 300)















