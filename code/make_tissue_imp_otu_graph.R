### Code for the tissue top 10 most important genera and OTUs Tax ID
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis", "stringr"))


# Read in Needed data
crc_unmatched_genera_occurances <- read_csv("data/process/tables/crc_RF_genera_unmatched_tissue_top10.csv")
crc_matched_genera_occurances <- read_csv("data/process/tables/adn_RF_genera_matched_tissue_top10.csv")
adn_genera_occurances <- read_csv("data/process/tables/adn_RF_genera_tissue_top10.csv")
crc_unmatched_otu_occurances <- read_csv("data/process/tables/crc_RF_otu_unmatched_tissue_top10.csv")
crc_matched_otu_occurances <- read_csv("data/process/tables/adn_RF_otu_matched_tissue_top10.csv")
adn_otu_occurances <- read_csv("data/process/tables/adn_RF_otu_tissue_top10.csv")


adn_tissue_sets <- 2
crc_unmatched_sets <- 4
crc_matched_sets <- 3

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################


adn_genera <- adn_genera_occurances %>% arrange(occurance) %>% 
  mutate(otu = str_replace(otu, "_", " "), 
         otu = factor(otu, 
                      levels = otu,
                      labels = otu), 
         occurance = occurance/adn_tissue_sets) %>% 
  ggplot(aes(otu, occurance)) + 
  geom_bar(stat = "identity", fill = '#8B8B83') + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, 1)) + 
  ggtitle("A") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.3, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




crc_unmatched_genera <- crc_unmatched_genera_occurances %>% arrange(occurance) %>% 
  mutate(otu = str_replace(otu, "_", " "), 
         otu = factor(otu, 
                      levels = otu,
                      labels = otu), 
         occurance = occurance/crc_unmatched_sets) %>% 
  ggplot(aes(otu, occurance)) + 
  geom_bar(stat = "identity", fill = "black") + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, 1)) + 
  ggtitle("B") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.35, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


crc_matched_genera <- crc_matched_genera_occurances %>% arrange(occurance) %>% 
  mutate(otu = str_replace(otu, "_", " "), 
         otu = factor(otu, 
                      levels = otu,
                      labels = otu), 
         occurance = occurance/crc_matched_sets) %>% 
  ggplot(aes(otu, occurance)) + 
  geom_bar(stat = "identity", fill = "black") + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, 1)) + 
  ggtitle("C") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.4, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())





adn_otu <- adn_otu_occurances %>% arrange(occurance) %>% 
  mutate(genus = str_replace(genus, "_", " "), 
         genus = factor(genus, 
                      levels = genus,
                      labels = genus), 
         occurance = occurance/adn_tissue_sets) %>% 
  ggplot(aes(genus, occurance)) + 
  geom_bar(stat = "identity", fill = '#8B8B83') + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, 1)) + 
  ggtitle("D") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.4, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())




crc_unmatched_otu <- crc_unmatched_otu_occurances %>% arrange(occurance) %>% 
  mutate(genus = str_replace_all(genus, "_", " "), 
         genus = factor(genus, 
                      levels = genus,
                      labels = genus), 
         occurance = occurance/crc_unmatched_sets) %>% 
  ggplot(aes(genus, occurance)) + 
  geom_bar(stat = "identity", fill = "black") + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, 1)) + 
  ggtitle("E") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.8, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


crc_matched_otu <- crc_matched_otu_occurances %>% arrange(occurance) %>% 
  mutate(genus = str_replace_all(genus, "_", " "), 
         genus = factor(genus, 
                      levels = genus,
                      labels = genus), 
         occurance = occurance/crc_matched_sets) %>% 
  ggplot(aes(genus, occurance)) + 
  geom_bar(stat = "identity", fill = "black") + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, 1)) + 
  ggtitle("F") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.8, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


tissue_graph <- grid.arrange(adn_genera, crc_unmatched_genera, crc_matched_genera, 
                            adn_otu, crc_unmatched_otu, crc_matched_otu, 
                            layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)))


ggsave("results/figures/FigureS7.pdf", 
       tissue_graph, width = 13, height = 8, dpi = 300)
