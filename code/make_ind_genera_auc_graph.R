### Code to create alpha RR graph Figure 2
### Stool specific measures only
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Load in needed data tables

stool_data <- read_csv("data/process/tables/ind_genera_auc_stool.csv") %>% 
  mutate(taxa = str_replace_all(taxa, "_unclassified", ""), 
         taxa = factor(taxa, 
                       levels = c("Clostridium_XI", "Ruminococcus", "Enterobacteriaceae", "Escherichia", "Fusobacterium", 
                                  "Parvimonas", "Peptostreptococcus", "Porphyromonas"), 
                       labels = c("Clostridium XI", "Ruminococcus", "Enterobacteriaceae", "Escherichia", "Fusobacterium", 
                                  "Parvimonas", "Peptostreptococcus", "Porphyromonas"))) %>% 
  
  mutate(auc = ifelse(taxa == "Ruminococcus" | taxa == "Clostridium XI", invisible(1-auc), invisible(auc)))

combined_stool <- stool_data %>% group_by(taxa) %>% 
  summarise(auc = median(auc)) %>% 
  mutate(study = "median") %>% ungroup() %>% 
  mutate(taxa = str_replace_all(taxa, "_unclassified", ""), 
         taxa = factor(taxa, 
                       levels = c("Clostridium XI", "Ruminococcus", "Enterobacteriaceae", "Escherichia", "Fusobacterium", 
                                  "Parvimonas", "Peptostreptococcus", "Porphyromonas"), 
                       labels = c("Clostridium XI", "Ruminococcus", "Enterobacteriaceae", "Escherichia", "Fusobacterium", 
                                  "Parvimonas", "Peptostreptococcus", "Porphyromonas")))

  

unmatched_tissue_data <- read_csv("data/process/tables/ind_genera_auc_unmatched_tissue.csv") %>% 
  mutate(taxa = str_replace_all(taxa, "_unclassified", ""), 
         taxa = factor(taxa, 
                       levels = c("Dorea", "Weissella", "Blautia"), 
                       labels = c("Dorea", "Weissella", "Blautia"))) %>% 
  mutate(auc = ifelse(taxa == "Dorea" | taxa == "Blautia", invisible(1-auc), invisible(auc)))

combined_unmatched_tissue <- unmatched_tissue_data %>% group_by(taxa) %>% 
  summarise(auc = median(auc)) %>% 
  mutate(study = "median") %>% ungroup() %>% 
  mutate(taxa = str_replace_all(taxa, "_unclassified", ""), 
          taxa = factor(taxa, 
                        levels = c("Dorea", "Weissella", "Blautia"), 
                        labels = c("Dorea", "Weissella", "Blautia")))

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################


stool_graph <- stool_data %>% 
  mutate(study = factor(study, 
                 levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                 labels = c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller"))) %>% 
  ggplot(aes(auc, taxa, color = study)) + 
  geom_vline(xintercept = 0.5, linetype = "dashed") + 
  geom_point(show.legend = T) + 
  geom_point(data = combined_stool, aes(auc, taxa), color = "black", size = 3.5, show.legend = F) + 
  coord_cartesian(xlim = c(0, 1)) + 
  theme_bw() + 
  scale_color_manual(name = "", values = c('#B0C4DE', '#8968CD', '#ED9121', '#006400', '#97D83FFF', '#8B4513', '#FDE725FF')) + 
  labs(x = "AUC", y = "") + ggtitle("A") + 
  annotate("text", label = paste("Carcinoma (Stool)"), x = 0.15, y = 8.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 10), 
        axis.text.x = element_text(size = 10))
  



tissue_graph <- unmatched_tissue_data %>% 
  mutate(study = factor(study, 
                        levels = c("burns", "chen",  "flemer", "sana"), 
                        labels = c("Burns", "Chen", "Flemer", "Sanapareddy"))) %>% 
  ggplot(aes(auc, taxa, color = study)) + 
  geom_vline(xintercept = 0.5, linetype = "dashed") + 
  geom_point(show.legend = T) + 
  geom_point(data = combined_unmatched_tissue, aes(auc, taxa), color = "black", size = 3.5, show.legend = F) + 
  theme_bw() + 
  scale_color_manual(name = "", values = c('#453581FF', '#CD6889', '#ED9121', '#8EE5EE'), 
                     guide = guide_legend(nrow = 2, ncol = 2)) + 
  labs(x = "AUC", y = "") + ggtitle("B") + 
  annotate("text", label = paste("Carcinoma\n(Unmatched Tissue)"), x = 0.2, y = 3.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 10), 
        axis.text.x = element_text(size = 10)) + 
  coord_fixed(xlim = c(0, 1), ratio = 0.3)


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

combined_graph <- grid.arrange(stool_graph, tissue_graph, ncol = 2, nrow = 1)

ggsave("results/figures/Figure3.pdf", combined_graph, width = 8.5, height = 4, dpi = 300)



