### Code to graph the estimated power for adenoma and carcinoma
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Load in the data
adn_power <- read_csv("data/process/tables/adn_predicted_pwr_and_n.csv")
crc_power <- read_csv("data/process/tables/cancer_predicted_pwr_and_n.csv")

combined_data <- adn_power %>% 
  filter(study %in% c("lu", "brim")) %>% 
  bind_rows(crc_power)

##############################################################################################
############################## List of functions to be used  #################################
##############################################################################################


### Study Colors by Viridis 
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


adn_study_power <- adn_power %>% 
  mutate(study = factor(study, 
                        levels = c("zeller", "lu", "hale", "flemer_t", "brim", "baxter"), 
                        labels = c("Zeller", "Lu", "Hale", "Flemer\n(Tissue)", "Brim", "Baxter")), 
         effect_size = factor(effect_size, 
                              levels = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), 
                              labels = c("1%", "5%", "10%", "15%", "20%", "25%", "30%"))) %>% 
  ggplot(aes(effect_size, study_power, color = study, group = study)) + 
  geom_jitter(width = 0.3, size = 3, show.legend = F) + 
  coord_cartesian(ylim = c(0,1)) + 
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), color = "gray") + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") + 
  labs(x = "Effect Size", y = "Study Power") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#FDE725FF', '#8B7500', '#006400', '#ED9121', 
                                '#34618DFF', '#8968CD')) + 
  annotate("text", label = paste("Adenoma"), x = 0.78, y = 1.03, size = 2.5) +
  theme(plot.title = element_text(face="bold", hjust = -0.15, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_study_power <- crc_power %>% 
  mutate(study = factor(study, 
                        levels = c("zeller", "weir", "wang", "sana", "hale", "geng", 
                                   "flemer", "dejea", "burns", "baxter", "ahn"), 
                        labels = c("Zeller", "Weir", "Wang", "Sanapareddy", "Hale", "Geng", 
                                   "Flemer","Dejea", "Burns", "Baxter", "Ahn")), 
         effect_size = factor(effect_size, 
                              levels = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), 
                              labels = c("1%", "5%", "10%", "15%", "20%", "25%", "30%"))) %>% 
  ggplot(aes(effect_size, study_power, color = study, group = study)) + 
  geom_jitter(width = 0.3, size = 3, show.legend = F) + 
  coord_cartesian(ylim = c(0,1)) + 
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), color = "gray") + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") + 
  labs(x = "Effect Size", y = "Study Power") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#FDE725FF', '#8B4513', '#97D83FFF', '#8EE5EE', 
                                '#006400', '#EEDC82', '#ED9121', '#1874CD', 
                                '#453581FF', '#8968CD', '#B0C4DE')) + 
  annotate("text", label = paste("Carcinoma"), x = 0.85, y = 1.03, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.15, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


combined_sample_n <- combined_data %>% 
  mutate(study = factor(study, 
                        levels = c("zeller", "weir", "wang", "sana", "lu", "hale", "geng", 
                                   "flemer", "dejea", "burns", "brim", "baxter", "ahn"), 
                        labels = c("Zeller", "Weir", "Wang", "Sanapareddy", "Lu", "Hale", "Geng", 
                                   "Flemer","Dejea", "Burns", "Brim", "Baxter", "Ahn")), 
         effect_size = factor(effect_size, 
                              levels = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), 
                              labels = c("1%", "5%", "10%", "15%", "20%", "25%", "30%"))) %>% 
  ggplot(aes(effect_size, log2(pwr80_needed_n), color = study, group = study)) + 
  geom_jitter(width = 0.3, size = 3) + 
  coord_cartesian(ylim = c(0,16)) + 
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), color = "gray") + 
  labs(x = "Effect Size", y = expression(Log["2"]~Sample~Number)) + theme_bw() + ggtitle("C") + 
  scale_color_manual(name = "Study", 
                     values = c('#FDE725FF', '#8B4513', '#97D83FFF', '#8EE5EE', '#8B7500',  
                                '#006400', '#EEDC82', '#ED9121', '#1874CD', '#453581FF', 
                                '#34618DFF', '#8968CD', '#B0C4DE')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.18, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############################## Execution of needed functions  ################################
##############################################################################################

power_graph <- grid.arrange(adn_study_power, crc_study_power, combined_sample_n, 
                            layout_matrix = rbind(c(1, 3), c(2, 3)))

ggsave("results/figures/Figure7.pdf", 
       power_graph, width = 10, height = 8, dpi = 300)




