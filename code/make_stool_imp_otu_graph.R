### Code for the stool top 10 most important genera and OTUs Tax ID
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis", "stringr"))

# Read in Needed data
stool_mda <- read_csv("data/process/tables/crc_RF_genera_stool_top10_mda.csv") %>% 
  mutate(otu = str_replace_all(otu, "_unclassified", ""), 
         otu = str_replace_all(otu, "\\.Shigella", ""), 
         otu = str_replace_all(otu, "_", " ")) %>% 
  group_by(study) %>% mutate(zscore = scale(mda_median)) %>% 
  ungroup()

stool_otu_mda <- read_csv("data/process/tables/crc_RF_otu_stool_top10_mda.csv") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "\\.Shigella", ""), 
         genus = str_replace_all(genus, "_", " ")) %>% 
  group_by(study, genus) %>% filter(mda_median == max(mda_median)) %>% 
  ungroup() %>% 
  group_by(study) %>% mutate(zscore = scale(mda_median)) %>% 
  ungroup()


crc_stool_sets <- 7


##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################


crc_genera <- stool_mda %>% 
  mutate(study = factor(study, 
                        levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                        labels= c(c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller")))) %>% 
  ggplot(aes(study, otu, fill = zscore)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient2(name = "Z-Score Median MDA", low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_bw() + ggtitle("A") + 
  labs(x = "", y = "") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 3.5), 
        plot.title = element_text(face="bold", hjust = -0.3, size = 20))


crc_otu <- stool_otu_mda %>% 
  mutate(study = factor(study, 
                        levels = c("ahn", "baxter", "flemer", "hale", "wang", "weir", "zeller"), 
                        labels= c(c("Ahn", "Baxter", "Flemer", "Hale", "Wang", "Weir", "Zeller")))) %>% 
  ggplot(aes(study, genus, fill = zscore)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient2(name = "Z-Score Median MDA", low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_bw() + ggtitle("B") + 
  labs(x = "", y = "") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        axis.text.y = element_text(face = "italic", size = 4.5), 
        plot.title = element_text(face="bold", hjust = -0.3, size = 20))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


stool_graph <- grid.arrange(crc_genera, crc_otu, nrow = 1, ncol = 2)


ggsave("results/figures/FigureS3.pdf", 
       stool_graph, width = 10, height = 11, dpi = 300)




