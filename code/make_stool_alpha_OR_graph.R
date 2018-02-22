### Code to create alpha RR graph Figure 2
### Stool specific measures only
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))


# Load in needed data tables
adn_all_stool <- read_csv("data/process/tables/alpha_adn_OR_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_adn_OR_ind_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V1-V3", 3), rep("V4", 6), rep("V3-V5", 3)))


crc_all_stool <- read_csv("data/process/tables/alpha_OR_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_OR_ind_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V3", 3), rep("V4", 3), rep("V3-V4", 3), 
                    rep("V4", 6), rep("V3-V5", 3), rep("V3-V4", 3)))


# Change V4 FFB90F and V3

### Used the `viridis_pal()(7)` to choose colors to use
#### zeller = V4 - #FFB90F
#### hale = V3-V5 - #443A83FF
#### brim = V1-V3 - #31688EFF
#### baxter = V4 - #FFB90F
#### weir = V4 - #FFB90F
#### wang = V3 - #FFC1C1
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

adn_stool_graph <- adn_all_stool %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "zeller", "hale", "brim", "baxter"), 
                        labels = c( "Pooled", "Zeller", "Hale", "Brim", "Baxter")), 
         region = factor(region, 
                         levels = c("combined", "V4", "V3-V5", "V1-V3"), 
                         labels = c("Combined", "V4", "V3-V5", "V1-V3")), 
         measure = factor(measure, 
                          levels = c("sobs", "shannoneven", "shannon"), 
                          labels = c("Observed OTUs", "Evenness", "Shannon Diversity"))) %>% 
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-5.2, 5.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = T) + 
  geom_point(size = 3, show.legend = T) + 
  facet_grid(. ~ measure) + 
  labs(x = expression(Log["2"]~Odds~Ratio), y = "") + theme_bw() + ggtitle("A") + 
  scale_color_manual(name = "Variable Region", values = c( '#000000', '#FFB90F', '#443A83FF', '#31688EFF')) + 
  annotate("text", label = paste("Adenoma (Stool)"), x = -2.75, y = 5.5, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))



crc_stool_graph <- crc_all_stool %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "zeller", "weir", "wang", "hale", 
                                   "flemer", "baxter", "ahn"), 
                        labels = c( "Pooled", "Zeller", "Weir", "Wang",  "Hale", 
                                    "Flemer", "Baxter", "Ahn")), 
         region = factor(region, 
                         levels = c("combined", "V3", "V4", "V3-V4", "V3-V5"), 
                         labels = c("Combined", "V3", "V4", "V3-V4", "V3-V5")), 
         measure = factor(measure, 
                          levels = c("sobs", "shannoneven", "shannon"), 
                          labels = c("Observed OTUs", "Evenness", "Shannon Diversity"))) %>%  
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-5.2, 5.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = T) + 
  geom_point(size = 3, show.legend = T) + 
  facet_grid(. ~ measure) + 
  labs(x = expression(Log["2"]~Odds~Ratio), y = "") + theme_bw() + ggtitle("B") + 
  scale_color_manual(name = "Variable Region", 
                     values = c('#000000', '#FFC1C1', '#FFB90F', '#35B779FF', '#443A83FF')) + 
  annotate("text", label = paste("Carcinoma (Stool)"), x = -2.70, y = 8.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

stool_alpha_RR <- grid.arrange(adn_stool_graph, crc_stool_graph)

ggsave("results/figures/Figure2.pdf", stool_alpha_RR, width = 8.5, height = 7, dpi = 300)

