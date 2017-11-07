### Code to measure AUCs of genus between studies within disease (Adenoma or Carcinoma)
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))

##############################################################################################
############################## List of functions to be used  #################################
##############################################################################################

# Function to generate the tables to be used
make_table <- function(studies, path_to_file, first_data_part_name, 
                       second_data_part_name, ending){
  
  tempFull <- sapply(studies, 
                     function(x) read_csv(paste(path_to_file, first_data_part_name, x, 
                                                ending, sep = "")) %>% 
                       mutate(model = x, model_type = "full"), simplify = F) %>% bind_rows()
  
  tempSelect <- sapply(studies, 
                       function(x) read_csv(paste(path_to_file, second_data_part_name, x, 
                                                  ending, sep = "")) %>% 
                         mutate(model = x, model_type = "select"), simplify = F) %>% bind_rows()
  
  tempALL <- tempFull %>% bind_rows(tempSelect)
  
  return(tempALL)
  
}


##############################################################################################
############################## Execution of needed functions  ################################
##############################################################################################

# Load needed data tables (adenoma -- stool)
adn_stool_studies <- c("baxter", "brim", "hale", "zeller")

adn_all_stool <- make_table(adn_stool_studies, "data/process/tables/", 
                            "adn_genus_stool_RF_full_", "adn_genus_stool_RF_select_", 
                            "_pvalue_summary.csv")

# Load needed data (adenoma -- tissue)
adn_tissue_studies <- c("flemer", "lu")

adn_all_tissue <- make_table(adn_tissue_studies, "data/process/tables/", 
                            "adn_genus_unmatched_tissue_RF_full_", 
                            "adn_genus_unmatched_tissue_RF_select_", 
                            "_pvalue_summary.csv")


# Load needed data (carcinoma -- stool)
crc_stool_studies <- c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller")

crc_all_stool <- make_table(crc_stool_studies, "data/process/tables/", 
                             "genus_stool_RF_full_", 
                             "genus_stool_RF_select_", 
                             "_pvalue_summary.csv")

# Load needed data (carcinoma -- tissue)
crc_matched_tissue_studies <- c("burns", "dejea", "geng")

crc_all_matched_tissue <- make_table(crc_matched_tissue_studies, "data/process/tables/", 
                            "genus_matched_tissue_RF_", 
                            "genus_matched_tissue_RF_select_", 
                            "_pvalue_summary.csv")


crc_unmatched_tissue_studies <- c("burns", "chen", "flemer", "sana")

crc_all_unmatched_tissue <- make_table(crc_unmatched_tissue_studies, "data/process/tables/", 
                                     "genus_unmatched_tissue_RF_", 
                                     "genus_unmatched_tissue_RF_select_", 
                                     "_pvalue_summary.csv")



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

adn_tissue_graph <- adn_all_tissue %>% 
  mutate(model = factor(model, 
                        levels = c("flemer", "lu"), 
                        labels = c("Flemer", "Lu")), 
         model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("flemer", "lu"), 
                        labels = c("Flemer", "Lu"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_point(size = 3.5, show.legend = F) + facet_grid(. ~ model_type) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + 
  scale_color_manual(name = "Study", 
                     values = c('#440154FF', '#FDE725FF')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


adn_stool_graph <- adn_all_stool %>% 
  mutate(model = factor(model, 
                        levels = c("baxter", "brim", "hale", "zeller"), 
                        labels = c("Baxter", "Brim", "Hale", "Zeller")), 
         model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("baxter", "brim", "hale", "zeller"), 
                        labels = c("Baxter", "Brim", "Hale", "Zeller"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_jitter(width = 0.2, size = 3.5, show.legend = F) + facet_grid(. ~ model_type) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#481D6FFF', '#34618DFF', '#67CC5CFF', '#FDE725FF')) + 
  annotate("text", label = paste("Adenoma (Stool)"), x = 1.0, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.09, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_unmatched_tissue_graph <- crc_all_unmatched_tissue %>% 
  mutate(model = factor(model, 
                        levels = c("burns", "chen", "flemer", "sana"), 
                        labels = c("Burns", "Chen", "Flemer", "Sanapareddy")), 
         model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("burns", "chen", "flemer", "sana"), 
                        labels = c("Burns", "Chen", "Flemer", "Sanapareddy"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_jitter(width = 0.25, size = 3.5, show.legend = F) + facet_grid(. ~ model_type) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#453581FF', '#3D4D8AFF', '#440154FF', '#1F998AFF')) + 
  annotate("text", label = paste("Carcinoma (Unmatched Tissue)"), x = 1.6, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_matched_tissue_graph <- crc_all_matched_tissue %>% 
  mutate(model = factor(model, 
                        levels = c("burns", "dejea", "geng"), 
                        labels = c("Burns", "Dejea", "Geng")), 
         model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("burns", "dejea", "geng"), 
                        labels = c("Burns", "Dejea", "Geng"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_jitter(width = 0.25, size = 3.5, show.legend = F) + facet_grid(. ~ model_type) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#453581FF', '#2B748EFF', '#CBE11EFF')) + 
  annotate("text", label = paste("Carcinoma (Matched Tissue)"), x = 1.3, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.text = element_text(size = 6), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_stool_graph <- crc_all_stool %>% 
  mutate(model = factor(model, 
                        levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                        labels = c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller")), 
         model_type = factor(model_type, 
                             levels = c("full", "select"), 
                             labels = c("All Genera", "CRC Associated\nGenera Only")), 
         study = factor(study, 
                        levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                        labels = c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_jitter(width = 0.2, size = 3.5, show.legend = F) + facet_grid(. ~ model_type) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#40BC72FF', '#481D6FFF', '#440154FF', '#67CC5CFF', 
                                '#97D83FFF', '#24878EFF', '#FDE725FF')) + 
  annotate("text", label = paste("Carcinoma (Stool)"), x = 1.4, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.09, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

crc_tissue_auc_graph <- grid.arrange(crc_unmatched_tissue_graph, crc_matched_tissue_graph)
stool_auc_graph <- grid.arrange(adn_stool_graph, crc_stool_graph)

ggsave("results/figures/FigureS5.pdf", 
       crc_tissue_auc_graph, width = 6, height = 6, dpi = 300)

ggsave("results/figures/FigureS4.pdf", 
       stool_auc_graph, width = 7, height = 6, dpi = 300)

ggsave("results/figures/FigureS6.pdf", 
       adn_tissue_graph, width = 3, height = 3, dpi = 300)

