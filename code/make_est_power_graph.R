### Code to graph the estimated power for adenoma and carcinoma
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Load in the data
adn_power <- read_csv("data/process/tables/adn_predicted_pwr_and_n.csv")
crc_power <- read_csv("data/process/tables/cancer_predicted_pwr_and_n.csv")



##############################################################################################
############################## List of functions to be used  #################################
##############################################################################################


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


adn_study_power <- adn_power %>% 
  mutate(study = factor(study, 
                        levels = c("zeller", "lu", "hale", "flemer_t", "brim", "baxter"), 
                        labels = c("Zeller", "Lu", "Hale", "Flemer\n(Tissue)", "Brim", "Baxter")), 
         effect_size = factor(effect_size, 
                              levels = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), 
                              labels = c("1%", "5%", "10%", "15%", "20%", "25%", "30%"))) %>% 
  ggplot(aes(effect_size, study_power, color = study, group = study)) + 
  geom_jitter(width = 0.3, size = 3) + 
  coord_cartesian(ylim = c(0,1)) + 
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), color = "gray") + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") + 
  labs(x = "Effect Size", y = "Study Power") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Study", 
                     values = c('#FDE725FF', '#FDE725FF', '#67CC5CFF', '#440154FF', 
                                '#34618DFF', '#481D6FFF')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


adn_sample_n <- adn_power %>% 
  mutate(study = factor(study, 
                        levels = c("zeller", "lu", "hale", "flemer_t", "brim", "baxter"), 
                        labels = c("Zeller", "Lu", "Hale", "Flemer\n(Tissue)", "Brim", "Baxter")), 
         effect_size = factor(effect_size, 
                              levels = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3), 
                              labels = c("1%", "5%", "10%", "15%", "20%", "25%", "30%"))) %>% 
  ggplot(aes(effect_size, log(pwr80_needed_n), color = study, group = study)) + 
  geom_jitter(width = 0.3, size = 3) + 
  coord_cartesian(ylim = c(0,12)) + 
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), color = "gray") + 
  labs(x = "Effect Size", y = expression(Log["10"]~Sample~Number)) + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Study", 
                     values = c('#FDE725FF', '#FDE725FF', '#67CC5CFF', '#440154FF', 
                                '#34618DFF', '#481D6FFF')) + 
  theme(plot.title = element_text(face="bold", hjust = -0.1, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))








##############################################################################################
############################## Execution of needed functions  ################################
##############################################################################################