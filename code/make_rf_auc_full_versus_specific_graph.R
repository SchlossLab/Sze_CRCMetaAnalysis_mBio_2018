### Code to measure AUCs of full data set and select genera dataset
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Load needed data tables (adenoma)
adn_tissue_matched <- 
  read_csv("data/process/tables/adn_genus_matched_tissue_RF_fullvsselect_pvalue_summary.csv") %>% 
  mutate(type = "matched")
adn_tissue_unmatched <- 
  read_csv("data/process/tables/adn_genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv") %>% 
  mutate(type = "unmatched")

adn_all_stool <- read_csv("data/process/tables/adn_genus_stool_RF_fullvsselect_pvalue_summary.csv") %>% 
  gather(key = model_type, value = AUC, full_model, select_model)
adn_all_tissue <- adn_tissue_unmatched %>% bind_rows(adn_tissue_matched) %>% 
  gather(key = model_type, value = AUC, full_model, select_model)

# Load in needed data tables (carcinoma)
crc_tissue_matched <- 
  read_csv("data/process/tables/genus_matched_tissue_RF_fullvsselect_pvalue_summary.csv") %>% 
  mutate(type = "matched")
crc_tissue_unmatched <- 
  read_csv("data/process/tables/genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv") %>% 
  mutate(type = "unmatched")

crc_all_stool <- read_csv("data/process/tables/genus_stool_RF_fullvsselect_pvalue_summary.csv") %>% 
  gather(key = model_type, value = AUC, full_model, select_model)
crc_all_tissue <- crc_tissue_unmatched %>% bind_rows(crc_tissue_matched) %>% 
  gather(key = model_type, value = AUC, full_model, select_model)

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

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################


adn_tissue_graph <- adn_all_tissue %>% 
  mutate(type = factor(type, 
                       levels = c("unmatched", "matched"), 
                       labels = c("Unmatched Tissue", "Matched Tissue")), 
         model_type = factor(model_type, 
                             levels = c("full_model", "select_model"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("flemer", "lu"), 
                        labels = c("Flemer", "Lu"))) %>% 
  ggplot(aes(model_type, AUC, color = study, group = type)) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) +
  geom_point(size = 3.5, show.legend = F) + facet_grid(. ~ type) + coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#440154FF', '#FDE725FF')) + 
  annotate("text", label = paste("Adenoma"), x = 0.6, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


adn_stool_graph <- adn_all_stool %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full_model", "select_model"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("baxter", "brim", "hale", "zeller"), 
                        labels = c("Baxter", "Brim", "Hale", "Zeller"))) %>% 
  ggplot(aes(model_type, AUC, color = study)) + 
  geom_jitter(width = 0.2, size = 3.5, show.legend = F) + coord_cartesian(ylim = c(0, 1.05)) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#481D6FFF', '#1F998AFF', '#67CC5CFF', '#FDE725FF')) + 
  annotate("text", label = paste("Adenoma (Stool)"), x = 0.80, y = 1.09, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.25, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))



crc_tissue_graph <- crc_all_tissue %>% 
  mutate(type = factor(type, 
                       levels = c("unmatched", "matched"), 
                       labels = c("Unmatched Tissue", "Matched Tissue")), 
         model_type = factor(model_type, 
                             levels = c("full_model", "select_model"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("burns", "chen", "dejea", "flemer", "geng", "sana"), 
                        labels = c("Burns", "Chen", "Dejea", "Flemer", "Geng", "Sanapareddy"))) %>% 
  ggplot(aes(model_type, AUC, color = study, group = type)) + 
  geom_jitter(width = 0.2, size = 3.5, show.legend = F) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  facet_grid(. ~ type) + coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#453581FF', '#3D4D8AFF', '#2B748EFF', '#440154FF', '#CBE11EFF', '#1F998AFF')) + 
  annotate("text", label = paste("Carcinoma"), x = 0.6, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_stool_graph <- crc_all_stool %>% 
  mutate(model_type = factor(model_type, 
                             levels = c("full_model", "select_model"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                        labels = c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller"))) %>% 
  ggplot(aes(model_type, AUC, color = study)) + 
  geom_jitter(width = 0.2, size = 3.5, show.legend = F) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#40BC72FF', '#481D6FFF', '#440154FF', '#67CC5CFF', '#97D83FFF', '#24878EFF', '#FDE725FF')) + 
  annotate("text", label = paste("Carcinoma (Stool)"), x = 0.8, y = 1.09, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.25, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

tissue_auc_graph <- grid.arrange(adn_tissue_graph, crc_tissue_graph)
stool_auc_graph <- grid.arrange(adn_stool_graph, crc_stool_graph, ncol = 2)

ggsave("results/figures/FigureS3.pdf", 
       tissue_auc_graph, width = 6, height = 6, dpi = 300)

ggsave("results/figures/FigureS2.pdf", 
       stool_auc_graph, width = 6, height = 6, dpi = 300)

