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
  geom_point(size = 3.5, show.legend = T) + facet_grid(. ~ type) + coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#ED9121', '#8B7500')) + 
  annotate("text", label = paste("Adenoma"), x = 0.6, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_blank(), 
        legend.background = element_rect(color = "black"), 
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
  geom_jitter(width = 0.2, size = 3.5, show.legend = T) + coord_cartesian(ylim = c(0, 1.05)) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#8968CD', '#34618DFF', '#006400', '#FDE725FF')) + 
  annotate("text", label = paste("Adenoma (Stool)"), x = 0.80, y = 1.09, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.25, size = 20), 
        legend.position = c(0.80, 0.22), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8), 
        legend.background = element_rect(color = "black"), 
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
  geom_jitter(width = 0.2, size = 3.5, show.legend = T) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  facet_grid(. ~ type) + coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#453581FF', '#CD6889', '#1874CD', '#ED9121', '#EEDC82', '#8EE5EE')) + 
  annotate("text", label = paste("Carcinoma"), x = 0.6, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_blank(), 
        legend.background = element_rect(color = "black"), 
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
  geom_jitter(width = 0.2, size = 3.5, show.legend = T) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "black", geom = "crossbar", size = 0.5, width = 0.5) + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "Model AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#B0C4DE', '#8968CD', '#ED9121', '#006400', '#97D83FFF', '#8B4513', '#FDE725FF')) + 
  annotate("text", label = paste("Carcinoma (Stool)"), x = 0.8, y = 1.09, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.25, size = 20), 
        legend.position = c(0.80, 0.22), 
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8), 
        legend.background = element_rect(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

tissue_auc_graph <- grid.arrange(adn_tissue_graph, crc_tissue_graph)
stool_auc_graph <- grid.arrange(adn_stool_graph, crc_stool_graph, ncol = 2)

ggsave("results/figures/FigureS3.pdf", 
       tissue_auc_graph, width = 6, height = 8, dpi = 300)

ggsave("results/figures/FigureS2.pdf", 
       stool_auc_graph, width = 6, height = 6, dpi = 300)

