### Code for the stool top 10 most important genera and OTUs Tax ID
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis", "stringr"))

# Read in Needed data
crc_genera_occurances <- read_csv("data/process/tables/crc_RF_genera_stool_top10.csv") %>% 
  mutate(otu = str_replace_all(otu, "_unclassified", ""), 
         otu = str_replace_all(otu, "\\.Shigella", ""))
crc_otu_occurances <- read_csv("data/process/tables/crc_RF_otu_stool_top10.csv") %>% 
  mutate(genus = str_replace_all(genus, "_unclassified", ""), 
         genus = str_replace_all(genus, "\\/Shigella", ""))

crc_stool_sets <- 6


##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################


crc_genera <- crc_genera_occurances %>% arrange(occurance) %>% 
  mutate(otu = str_replace(otu, "_", " "), 
         otu = factor(otu, 
                      levels = otu,
                      labels = otu), 
         occurance = occurance) %>% 
  ggplot(aes(otu, occurance)) + 
  geom_bar(stat = "identity", fill = '#DB7093') + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, 6)) + 
  ggtitle("A") + 
  scale_y_continuous(expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.3, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


crc_otu <- crc_otu_occurances %>% arrange(occurance) %>% 
  mutate(genus = str_replace(genus, "_", " "), 
         genus = factor(genus, 
                        levels = genus,
                        labels = genus), 
         occurance = occurance) %>% 
  ggplot(aes(genus, occurance)) + 
  geom_bar(stat = "identity", fill = '#B0171F') + 
  theme_bw() + 
  labs(x = "", y = "Occurrence Across Studies") + 
  coord_flip(ylim = c(0, 6)) + 
  ggtitle("B") + 
  scale_y_continuous(expand = c(0,0.03)) + 
  theme(axis.text.y = element_text(face = "italic"), 
        plot.title = element_text(face="bold", hjust = -0.3, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


stool_graph <- grid.arrange(crc_genera, crc_otu, nrow = 1, ncol = 2)


ggsave("results/figures/Figure6.pdf", 
       stool_graph, width = 8, height = 6, dpi = 300)




