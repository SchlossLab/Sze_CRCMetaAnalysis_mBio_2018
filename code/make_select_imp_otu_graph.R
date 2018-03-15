### Code for the crc most imp taxa in select models
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis", "stringr"))

# Read in Needed data
stool_mda <- read_csv("data/process/tables/ALL_genus_stool_RF_select_imp_vars.csv") %>% 
  mutate(otu = str_replace_all(otu, "_unclassified", ""), 
         otu = str_replace_all(otu, "\\.Shigella", ""), 
         otu = str_replace_all(otu, "_", " ")) %>% 
  group_by(study) %>% mutate(zscore = scale(mda_median)) %>% 
  ungroup()

unmatched_tissue_mda <- read_csv("data/process/tables/ALL_genus_unmatched_tissue_RF_select_imp_vars.csv") %>% 
  mutate(otu = str_replace_all(otu, "_unclassified", ""), 
         otu = str_replace_all(otu, "\\.Shigella", ""), 
         otu = str_replace_all(otu, "_", " ")) %>% 
  group_by(study, otu) %>% filter(mda_median == max(mda_median)) %>% 
  ungroup() %>% 
  group_by(study) %>% mutate(zscore = scale(mda_median)) %>% 
  ungroup()


stool_vars <- stool_mda %>% group_by(otu) %>% summarise(median_zscore = mean(zscore)) %>% 
  arrange(median_zscore) %>% ungroup()

unmatched_tissue_vars <- unmatched_tissue_mda %>% group_by(otu) %>% summarise(median_zscore = mean(zscore)) %>% 
  arrange(median_zscore) %>% ungroup()


##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################


imp_stool_model <- stool_mda %>% 
  mutate(study = factor(study, 
                        levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                        labels= c(c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller"))), 
         otu = factor(otu, 
                      levels = stool_vars$otu, 
                      labels = stool_vars$otu)) %>% 
  ggplot(aes(study, otu, fill = zscore)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient2(name = "Z-Score MDA", low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_bw() + ggtitle("A") + 
  labs(x = "", y = "") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 10), 
        plot.title = element_text(face="bold", hjust = -0.3, size = 20))


imp_unmatched_tissue_model <- unmatched_tissue_mda %>% 
  mutate(study = factor(study, 
                        levels = c("burns", "chen", "flemer", "sana"), 
                        labels= c(c("Burns", "Chen", "Flemer", "Sanapareddy"))), 
         otu = factor(otu, 
                      levels = unmatched_tissue_vars$otu, 
                      labels = unmatched_tissue_vars$otu)) %>% 
  ggplot(aes(study, otu, fill = zscore)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient2(name = "Z-Score MDA", low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_bw() + ggtitle("B") + 
  labs(x = "", y = "") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 10), 
        plot.title = element_text(face="bold", hjust = -0.15, size = 20))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


stool_graph <- grid.arrange(imp_stool_model, imp_unmatched_tissue_model, nrow = 1, ncol = 2)


ggsave("results/figures/Figure4.pdf", 
       stool_graph, width = 11, height = 6, dpi = 300)




