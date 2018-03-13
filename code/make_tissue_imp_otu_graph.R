### Code for the tissue top 10 most important genera and OTUs Tax ID
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis", "stringr"))


# Read in Needed data
crc_unmatched_genera <- read_csv("data/process/tables/crc_RF_genera_unmatched_tissue_top10_mda.csv") %>% 
  mutate(otu = str_replace_all(otu, "_unclassified", ""), 
         otu = str_replace_all(otu, "\\.", "/"), 
         otu = str_replace_all(otu, "_", " "), 
         otu = str_replace_all(otu, "Gp", "Acidobacteria Gp"), 
         otu = str_replace_all(otu, "Acidobacteria Acidobacteria Gp", "Acidobacteria Gp")) %>% 
  group_by(study) %>% mutate(zscore = scale(mda_median)) %>% 
  ungroup()

crc_matched_genera <- read_csv("data/process/tables/adn_RF_genera_matched_tissue_top10_mda.csv") %>% 
  mutate(otu = str_replace_all(otu, "_unclassified", ""), 
         otu = str_replace_all(otu, "\\.", "/"), 
         otu = str_replace_all(otu, "_", " "), 
         otu = str_replace_all(otu, "Gp", "Acidobacteria Gp"), 
         otu = str_replace_all(otu, "Acidobacteria Acidobacteria Gp", "Acidobacteria Gp")) %>% 
  group_by(study) %>% mutate(zscore = scale(mda_median)) %>% 
  ungroup() %>% filter(otu != "NA/")


crc_unmatched_otu <- read_csv("data/process/tables/crc_RF_otu_unmatched_tissue_top10_mda.csv") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "\\.", "/"), 
         genus = str_replace_all(genus, "_", " "), 
         genus = str_replace_all(genus, "Gp", "Acidobacteria Gp"), 
         genus = str_replace_all(genus, "Acidobacteria Acidobacteria Gp", "Acidobacteria Gp")) %>% 
  group_by(study, genus) %>% filter(mda_median == max(mda_median)) %>% 
  ungroup() %>% 
  group_by(study) %>% mutate(zscore = scale(mda_median)) %>% 
  ungroup()

crc_matched_otu <- read_csv("data/process/tables/adn_RF_otu_matched_tissue_top10_mda.csv") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "\\.", "/"), 
         genus = str_replace_all(genus, "_", " "), 
         genus = str_replace_all(genus, "Gp", "Acidobacteria Gp"), 
         genus = str_replace_all(genus, "Acidobacteria Acidobacteria Gp", "Acidobacteria Gp")) %>% 
  group_by(study, genus) %>% filter(mda_median == max(mda_median)) %>% 
  ungroup() %>% 
  group_by(study) %>% mutate(zscore = scale(mda_median)) %>% 
  ungroup()


crc_unmatched_sets <- 4
crc_matched_sets <- 3

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

crc_matched_genera_graph <- crc_matched_genera %>% 
  mutate(study = factor(study, 
                        levels = c("burns", "dejea", "geng"), 
                        labels= c("Burns", "Dejea", "Geng"))) %>% 
  ggplot(aes(study, otu, fill = zscore)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient2(name = "Z-Score Median MDA", low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_bw() + ggtitle("A") + 
  labs(x = "", y = "") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 1.5), 
        plot.title = element_text(face="bold", hjust = -0.1, size = 20))


crc_unmatched_genera_graph <- crc_unmatched_genera %>% 
  mutate(study = factor(study, 
                        levels = c("burns", "chen", "flemer", "sana"), 
                        labels= c("Burns", "Chen", "Flemer", "Sanapareddy"))) %>% 
  ggplot(aes(study, otu, fill = zscore)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient2(name = "Z-Score Median MDA", low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_bw() + ggtitle("B") + 
  labs(x = "", y = "") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 1), 
        plot.title = element_text(face="bold", hjust = -0.1, size = 20))


crc_matched_otu_graph <- crc_matched_otu %>% 
  mutate(study = factor(study, 
                        levels = c("burns", "dejea", "geng"), 
                        labels= c("Burns", "Dejea", "Geng"))) %>% 
  ggplot(aes(study, genus, fill = zscore)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient2(name = "Z-Score Median MDA", low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_bw() + ggtitle("C") + 
  labs(x = "", y = "") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 3.5), 
        plot.title = element_text(face="bold", hjust = -0.25, size = 20))


crc_unmatched_otu_graph <- crc_unmatched_otu %>% 
  mutate(study = factor(study, 
                        levels = c("burns", "chen", "flemer", "sana"), 
                        labels= c("Burns", "Chen", "Flemer", "Sanapareddy"))) %>% 
  ggplot(aes(study, genus, fill = zscore)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient2(name = "Z-Score Median MDA", low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_bw() + ggtitle("D") + 
  labs(x = "", y = "") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 3.5), 
        plot.title = element_text(face="bold", hjust = -0.25, size = 20))



##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


tissue_graph <- grid.arrange(crc_matched_genera_graph, crc_unmatched_genera_graph, 
                             crc_matched_otu_graph, crc_unmatched_otu_graph,  
                            layout_matrix = rbind(c(1, 2), c(3, 4)))


ggsave("results/figures/FigureS4.pdf", 
       tissue_graph, width = 10, height = 13, dpi = 300)
