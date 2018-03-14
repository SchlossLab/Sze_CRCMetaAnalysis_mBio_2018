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
                       second_data_part_name, ending, select_included = T){
  
  
  tempFull <- sapply(studies, 
                     function(x) read_csv(paste(path_to_file, first_data_part_name, x, 
                                                ending, sep = "")) %>% 
                       mutate(model = x, model_type = "full"), simplify = F) %>% bind_rows()
  
  if(select_included == T){
    
    tempSelect <- sapply(studies, 
                         function(x) read_csv(paste(path_to_file, second_data_part_name, x, 
                                                    ending, sep = "")) %>% 
                           mutate(model = x, model_type = "select"), simplify = F) %>% bind_rows()
    
    tempALL <- tempFull %>% bind_rows(tempSelect)
  } else{
    
    tempALL <- tempFull
  }
  
  
  
  
  
  return(tempALL)
  
}


##############################################################################################
############################## Execution of needed functions  ################################
##############################################################################################

# Load needed data tables (adenoma -- stool)
adn_stool_studies <- c("baxter", "brim", "hale", "zeller")

adn_all_stool <- make_table(adn_stool_studies, "data/process/tables/", 
                            "adn_ALL_genus_stool_RF_full_", "adn_ALL_genus_stool_RF_full_", 
                            "_pvalue_summary.csv",select_included = F)

# Load needed data (adenoma -- tissue)
adn_tissue_studies <- c("flemer", "lu")

adn_all_tissue <- make_table(adn_tissue_studies, "data/process/tables/", 
                            "adn_ALL_genus_unmatched_tissue_RF_full_", 
                            "adn_ALL_genus_unmatched_tissue_RF_full_", 
                            "_pvalue_summary.csv", select_included = F)


# Load needed data (carcinoma -- stool)
crc_stool_studies <- c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller")

crc_all_stool <- make_table(crc_stool_studies, "data/process/tables/", 
                             "ALL_genus_stool_RF_full_", 
                             "ALL_genus_stool_RF_select_", 
                             "_pvalue_summary.csv")

# Load needed data (carcinoma -- tissue)
crc_matched_tissue_studies <- c("burns", "dejea", "geng")

crc_all_matched_tissue <- make_table(crc_matched_tissue_studies, "data/process/tables/", 
                            "ALL_genus_matched_tissue_RF_", 
                            "ALL_genus_matched_tissue_RF_", 
                            "_pvalue_summary.csv", select_included = F)


crc_unmatched_tissue_studies <- c("burns", "chen", "flemer", "sana")

crc_all_unmatched_tissue <- make_table(crc_unmatched_tissue_studies, "data/process/tables/", 
                                     "ALL_genus_unmatched_tissue_RF_", 
                                     "ALL_genus_unmatched_tissue_RF_select_", 
                                     "_pvalue_summary.csv")



##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

## Study Colors by Viridis 
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

adn_tissue_graph <- adn_all_tissue %>% 
  filter(study != model) %>% 
  mutate(model = factor(model, 
                        levels = c("flemer", "lu"), 
                        labels = c("Flemer", "Lu\n(Matched)")), 
         model_type = factor(model_type, 
                             levels = c("select", "full"), 
                             labels = c("Significant\nOR Taxa", "All Taxa")), 
         study = factor(study, 
                        levels = c("flemer", "lu"), 
                        labels = c("Flemer", "Lu\n(Matched)"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_point(size = 3.5, position = position_dodge(width = 0.3), show.legend = T) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "red", geom = "crossbar", size = 0.5, width = 0.5) + 
  facet_grid(. ~ model_type) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#ED9121', '#8B7500')) + 
  annotate("text", label = paste("Adenoma Tissue"), x = 0.8, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_blank(), 
        legend.background = element_rect(color = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


adn_stool_graph <- adn_all_stool %>% 
  filter(study != model) %>% 
  mutate(model = factor(model, 
                        levels = c("baxter", "brim", "hale", "zeller"), 
                        labels = c("Baxter", "Brim", "Hale", "Zeller")), 
         model_type = factor(model_type, 
                             levels = c("select", "full"), 
                             labels = c("Significant\nOR Taxa", "All Taxa")), 
         study = factor(study, 
                        levels = c("baxter", "brim", "hale", "zeller"), 
                        labels = c("Baxter", "Brim", "Hale", "Zeller"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_point(position = position_dodge(width = 0.5), size = 3.5, show.legend = T) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "red", geom = "crossbar", size = 0.5, width = 0.5) + 
  facet_grid(. ~ model_type) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#8968CD', '#34618DFF', '#006400', '#FDE725FF')) + 
  annotate("text", label = paste("Adenoma (Feces)"), x = 1.0, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.09, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_blank(), 
        legend.background = element_rect(color = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_unmatched_tissue_graph <- crc_all_unmatched_tissue %>% 
  filter(study != model) %>% 
  mutate(model = factor(model, 
                        levels = c("burns", "chen", "flemer", "sana"), 
                        labels = c("Burns", "Chen", "Flemer", "Sanapareddy")), 
         model_type = factor(model_type, 
                             levels = c("select", "full"), 
                             labels = c("Significant\nOR Taxa", "All Taxa")), 
         study = factor(study, 
                        levels = c("burns", "chen", "flemer", "sana"), 
                        labels = c("Burns", "Chen", "Flemer", "Sanapareddy"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_point(position = position_dodge(width = 0.3), size = 3.5, show.legend = T) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "red", geom = "crossbar", size = 0.5, width = 0.5) + 
  facet_grid(. ~ model_type) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + ggtitle("C") + 
  scale_color_manual(name = "Study", 
                     values = c('#453581FF', '#CD6889', '#ED9121', '#8EE5EE')) + 
  annotate("text", label = paste("Carcinoma (Unmatched Tissue)"), x = 1.6, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_blank(), 
        legend.background = element_rect(color = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_matched_tissue_graph <- crc_all_matched_tissue %>% 
  filter(study != model) %>% 
  mutate(model = factor(model, 
                        levels = c("burns", "dejea", "geng"), 
                        labels = c("Burns", "Dejea", "Geng")), 
         model_type = factor(model_type, 
                             levels = c("select", "full"), 
                             labels = c("Significant\nOR Taxa", "All Taxa")), 
         study = factor(study, 
                        levels = c("burns", "dejea", "geng"), 
                        labels = c("Burns", "Dejea", "Geng"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_point(position = position_dodge(width = 0.3), size = 3.5, show.legend = T) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "red", geom = "crossbar", size = 0.5, width = 0.5)  + 
  facet_grid(. ~ model_type) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#453581FF', '#1874CD', '#EEDC82')) + 
  annotate("text", label = paste("Carcinoma (Matched Tissue)"), x = 1.3, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_blank(), 
        legend.background = element_rect(color = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_stool_graph <- crc_all_stool %>% 
  filter(study != model) %>% 
  mutate(model = factor(model, 
                        levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                        labels = c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller")), 
         model_type = factor(model_type, 
                             levels = c("select", "full"), 
                             labels = c("Significant\nOR Taxa", "All Taxa")), 
         study = factor(study, 
                        levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                        labels = c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller"))) %>% 
  ggplot(aes(model, auc, color = study, group = model_type)) + 
  geom_point(position = position_dodge(width = 0.3), size = 3.5, show.legend = T) + 
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               colour = "red", geom = "crossbar", size = 0.5, width = 0.5) + 
  facet_grid(. ~ model_type) +
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  coord_cartesian(ylim = c(0, 1.05)) + 
  labs(x = "", y = "AUC") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#B0C4DE', '#8968CD', '#ED9121', '#006400', 
                                '#97D83FFF', '#8B4513', '#FDE725FF')) + 
  annotate("text", label = paste("Carcinoma (Feces)"), x = 1.4, y = 1.07, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.09, size = 20), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6),
        legend.title = element_blank(), 
        legend.background = element_rect(color = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

crc_tissue_auc_graph <- grid.arrange(adn_tissue_graph, crc_matched_tissue_graph, crc_unmatched_tissue_graph, 
                                     layout_matrix = rbind(c(1, 2), c(3, 3)))
stool_auc_graph <- grid.arrange(adn_stool_graph, crc_stool_graph)

ggsave("results/figures/FigureS5.pdf", 
       crc_tissue_auc_graph, width = 6, height = 8, dpi = 300)

ggsave("results/figures/Figure6.pdf", 
       stool_auc_graph, width = 7, height = 8, dpi = 300)


