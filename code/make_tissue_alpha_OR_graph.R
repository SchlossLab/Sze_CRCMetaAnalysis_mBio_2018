### Code to create alpha RR graph Figure 2
### Tissue specific measures only
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))

# Load in needed data tables
adn_all_tissue <- read_csv("data/process/tables/alpha_adn_OR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_adn_OR_ind_tissue_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V3-V4", 6)))


crc_all_tissue <- read_csv("data/process/tables/alpha_OR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_OR_ind_tissue_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V3-V5", 3), rep("V1-V2", 6), rep("V5-V6", 3), 
                    rep("V3-V4", 3), rep("V1-V3", 3)))



### Used the `viridis_pal()(7)` to choose colors to use
#### zeller = V4 - #440154FF
#### hale = V3-V5 - #443A83FF
#### brim = V1-V3 - #31688EFF
#### baxter = V4 - #440154FF
#### weir = V4 - #440154FF
#### wang = V3 - #21908CFF
### flemer = V3-V4 - #35B779FF
### ahn = V3-V4 - #35B779FF
### lu = V3-V4 - #35B779FF
### dejea = V3-V5 - #443A83FF
### geng = V1-V2 - #8FD744FF
### sana = V1-V2 - #8FD744FF
### burns = V5-V6 - #FDE725FF
### chen = V1-V3 - #31688EFF

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################


adn_tissue_graph <- adn_all_tissue %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "flemer", "lu"), 
                        labels = c( "Pooled", "Flemer", "Lu")), 
         region = factor(region, 
                         levels = c("combined", "V3-V4"), 
                         labels = c("Combined", "V3-V4")), 
         measure = factor(measure, 
                          levels = c("sobs", "shannoneven", "shannon"), 
                          labels = c("Observed OTUs", "Evenness", "Shannon Diversity")))  %>% 
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-6.2, 6.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = T) + 
  geom_point(size = 3, show.legend = T) + 
  facet_grid(. ~ measure) + 
  labs(x = expression(Log["2"]~Odds~Ratio), y = "") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Variable Region", values = c('#000000', '#35B779FF')) + 
  annotate("text", label = paste("Adenoma\n(Tissue)"), x = -4.70, y = 3.3, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_tissue_graph <- crc_all_tissue %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "sana", "geng", "flemer", "dejea", "chen", "burns"), 
                        labels = c( "Pooled", "Sanapareddy", "Geng", "Flemer", "Dejea", "Chen", "Burns")), 
         region = factor(region, 
                         levels = c("combined", "V1-V2", "V1-V3", "V3-V4", "V3-V5", "V5-V6"), 
                         labels = c("Combined", "V1-V2", "V1-V3", "V3-V4", "V3-V5", "V5-V6")), 
         measure = factor(measure, 
                          levels = c("sobs", "shannoneven", "shannon"), 
                          labels = c("Observed OTUs", "Evenness", "Shannon Diversity"))) %>%  
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-6.2, 6.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = T) + 
  geom_point(size = 3, show.legend = T) + 
  facet_grid(. ~ measure) + 
  labs(x = expression(Log["2"]~Odds~Ratio), y = "") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Variable Region", values = c('#000000', '#8FD744FF', '#31688EFF', 
                                '#35B779FF', '#443A83FF', '#FDE725FF')) + 
  annotate("text", label = paste("Carcinoma\n(Tissue)"), x = -4.4, y = 7, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.14, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

tissue_alpha_RR <- grid.arrange(adn_tissue_graph, crc_tissue_graph)

ggsave("results/figures/FigureS1.pdf", tissue_alpha_RR, width = 8.5, height = 7, dpi = 300)


