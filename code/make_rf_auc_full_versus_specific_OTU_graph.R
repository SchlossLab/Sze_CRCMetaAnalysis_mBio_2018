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
# flemer - #ED9121
# lu - #8B7500
# burns - #453581FF
# chen - #CD6889
# sana - #8EE5EE
# dejea - #1874CD
# geng - #EEDC82
# brim - #34618DFF
# zeller - #FDE725FF
# baxter - #8968CD
# hale - #006400
# wang - #97D83FFF
# weir - #8B4513
# ahn - #B0C4DE

# Creates a median secondary table for adenoma tissue
adn_tissue_medians <- adn_tissue %>% 
  filter(model == "act") %>% group_by(model_type) %>% 
  summarise(auc_median = median(AUC)) %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "Select Genera Only")))

# Creates the adenoma tissue graph
adn_tissue_graph <- adn_tissue %>% 
  mutate(type = factor(type, 
                       levels = c("unmatched", "matched"), 
                       labels = c("Unmatched Tissue", "Matched Tissue")), 
         model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "Select Genera Only")), 
         study = factor(study, 
                        levels = c("flemer", "lu"), 
                        labels = c("Flemer", "Lu\n(matched)"))) %>% 
  filter(grepl("rand", model) != T) %>% 
  ggplot(aes(study, AUC, color = study, group = study)) + 
  geom_hline(data = adn_tissue_medians, aes(yintercept = auc_median), linetype = "solid", color = "red", alpha = 0.7) + 
  geom_point(size = 3.5, show.legend = F) + 
  geom_errorbar(aes(ymin = AUC-deviation, ymax = AUC+deviation), 
                width = 0.25, size = 0.7, show.legend = F) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  facet_grid(. ~ model_type) + coord_flip(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#ED9121', '#8B7500')) + 
  annotate("text", label = paste("Adenoma"), x = 2.5, y = 0.05, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))

# Creates a median secondary table for adenoma stool
adn_stool_medians <- adn_all_stool %>% 
  filter(model == "act") %>% group_by(model_type) %>% 
  summarise(auc_median = median(AUC)) %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "Select Genera Only")))

# Creates the adenoma stool graph
adn_stool_graph <- adn_all_stool %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "Select Genera Only")), 
         study = factor(study, 
                        levels = c("baxter", "brim", "hale", "zeller"), 
                        labels = c("Baxter", "Brim", "Hale", "Zeller"))) %>% 
  filter(grepl("rand", model) != T) %>% 
  ggplot(aes(study, AUC, color = study, group = study)) + 
  geom_hline(data = adn_stool_medians, aes(yintercept = auc_median), linetype = "solid", color = "red", alpha = 0.7) + 
  geom_point(size = 3.5, show.legend = F) + 
  geom_errorbar(aes(ymin = AUC-deviation, ymax = AUC+deviation), 
                width = 0.25, size = 0.7, show.legend = F) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
   coord_flip(ylim = c(0, 1.05)) + 
  facet_grid(. ~ model_type) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#8968CD', '#1F998AFF', '#006400', '#FDE725FF')) + 
  annotate("text", label = paste("Adenoma (Stool)"), x = 4.5, y = 0.1, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


# Creates a median secondary table for carcinoma tissue
crc_tissue_medians <- crc_all_tissue %>% 
  filter(model == "act") %>% group_by(model_type) %>% 
  summarise(auc_median = median(AUC)) %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "Select Genera Only")))

# Creates the crc tissue graph
crc_tissue_graph <- crc_all_tissue %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "Select Genera Only")), 
         combined_study = factor(combined_study, 
                        levels = c("sana_unmatched", "geng_matched", "flemer_unmatched", 
                                   "dejea_matched", "chen_unmatched", "burns_unmatched",
                                   "burns_matched"), 
                        labels = c("Sanapareddy", "Geng\n(matched)", "Flemer", "Dejea\n(matched)", 
                                   "Chen", "Burns\n(matched)", "Burns"))) %>% 
  filter(grepl("rand", model) != T) %>% 
  ggplot(aes(combined_study, AUC, color = combined_study, group = combined_study)) + 
  geom_hline(data = crc_tissue_medians, aes(yintercept = auc_median), linetype = "solid", color = "red", alpha = 0.7) +
  geom_point(size = 3.5, show.legend = F) + 
  geom_errorbar(aes(ymin = AUC-deviation, ymax = AUC+deviation), 
                width = 0.25, size = 0.7, show.legend = F) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  facet_grid(. ~ model_type) + coord_flip(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#8EE5EE', '#EEDC82', '#ED9121', '#1874CD', 
                                '#CD6889', '#453581FF', '#453581FF')) + 
  annotate("text", label = paste("Carcinoma"), x = 7.35, y = 0.06, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.15, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


# Creates a median secondary table for carcinoma stool
crc_stool_medians <- crc_all_stool %>% 
  filter(model == "act") %>% group_by(model_type) %>% 
  summarise(auc_median = median(AUC)) %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "Select Genera Only")))

# Creates the crc stool graph
crc_stool_graph <- crc_all_stool %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("Full Community", "Select Genera Only")), 
         study = factor(study, 
                        levels = c("zeller", "weir", "wang", "hale", "flemer", "baxter", "ahn"), 
                        labels = c("Zeller", "Weir", "Wang", "Hale", "Flemer", "Baxter", "Ahn"))) %>% 
  filter(grepl("rand", model) != T) %>% 
  ggplot(aes(study, AUC, color = study, group = study)) + 
  geom_hline(data = crc_stool_medians, aes(yintercept = auc_median), linetype = "solid", color = "red", alpha = 0.7) +
  geom_point(size = 3.5, show.legend = F) + 
  geom_errorbar(aes(ymin = AUC-deviation, ymax = AUC+deviation), 
                width = 0.25, size = 0.7, show.legend = F) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_flip(ylim = c(0, 1.05)) + 
  facet_grid(. ~ model_type) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#FDE725FF', '#8B4513', '#97D83FFF', '#006400', 
                                '#ED9121', '#8968CD', '#B0C4DE')) + 
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















