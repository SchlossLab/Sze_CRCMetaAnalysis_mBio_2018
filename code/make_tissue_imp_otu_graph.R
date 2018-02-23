### Code for the tissue top 10 most important genera and OTUs Tax ID
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis", "stringr"))


# Read in Needed data
crc_unmatched_genera_occurances <- read_csv("data/process/tables/crc_RF_genera_unmatched_tissue_top10.csv") %>% 
  mutate(otu = str_replace_all(otu, "_unclassified", ""), 
         otu = str_replace_all(otu, "\\.", "/"))
crc_matched_genera_occurances <- read_csv("data/process/tables/adn_RF_genera_matched_tissue_top10.csv") %>% 
  mutate(otu = str_replace_all(otu, "_unclassified", ""), 
         otu = str_replace_all(otu, "\\.", "/"))
crc_unmatched_otu_occurances <- read_csv("data/process/tables/crc_RF_otu_unmatched_tissue_top10.csv") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "\\.", "/"))
crc_matched_otu_occurances <- read_csv("data/process/tables/adn_RF_otu_matched_tissue_top10.csv") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "\\.", "/"))


crc_unmatched_sets <- 4
crc_matched_sets <- 3

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################


crc_unmatched_genera <- crc_unmatched_genera_occurances %>% arrange(occurance) %>% 
  mutate(otu = str_replace(otu, "_", " "), 
         otu = factor(otu, 
                      levels = otu,
                      labels = otu), 
         occurance = occurance) %>% 
  ggplot(aes(otu, occurance)) + 
  geom_bar(stat = "identity", fill = '#DB7093') + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, crc_unmatched_sets)) + 
  ggtitle("B") + 
  scale_y_continuous(expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.35, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


crc_matched_genera <- crc_matched_genera_occurances %>% arrange(occurance) %>% 
  mutate(otu = str_replace(otu, "_", " "), 
         otu = factor(otu, 
                      levels = otu,
                      labels = otu), 
         occurance = occurance) %>% 
  ggplot(aes(otu, occurance)) + 
  geom_bar(stat = "identity", fill = '#DB7093') + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, crc_matched_sets)) + 
  ggtitle("A") + 
  scale_y_continuous(expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.20, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


crc_unmatched_otu <- crc_unmatched_otu_occurances %>% arrange(occurance) %>% 
  mutate(genus = str_replace_all(genus, "_", " "), 
         genus = factor(genus, 
                      levels = genus,
                      labels = genus), 
         occurance = occurance) %>% 
  ggplot(aes(genus, occurance)) + 
  geom_bar(stat = "identity", fill = '#B0171F') + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, crc_unmatched_sets)) + 
  ggtitle("D") + 
  scale_y_continuous(expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.35, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


crc_matched_otu <- crc_matched_otu_occurances %>% arrange(occurance) %>% 
  mutate(genus = str_replace_all(genus, "_", " "), 
         genus = factor(genus, 
                      levels = genus,
                      labels = genus), 
         occurance = occurance) %>% 
  ggplot(aes(genus, occurance)) + 
  geom_bar(stat = "identity", fill = '#B0171F') + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, crc_matched_sets)) + 
  ggtitle("C") + 
  scale_y_continuous(expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.35, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


tissue_graph <- grid.arrange(crc_matched_genera, crc_unmatched_genera, 
                             crc_matched_otu, crc_unmatched_otu,  
                            layout_matrix = rbind(c(1, 2), c(3, 4)))


ggsave("results/figures/FigureS4.pdf", 
       tissue_graph, width = 13, height = 8, dpi = 300)
