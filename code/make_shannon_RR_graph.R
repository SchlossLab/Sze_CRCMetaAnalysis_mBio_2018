### Code to create needed alpha diversity RR graphs
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra"))

# Load needed data tables (adenoma)
adn_all_stool <- read_csv("data/process/tables/alpha_adn_RR_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_adn_RR_ind_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V1-V3", 3), rep("V4", 6), rep("V3-V5", 3)))

adn_all_tissue <- read_csv("data/process/tables/alpha_adn_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_adn_RR_ind_tissue_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V3-V4", 6)))

# Load in needed data tables (carcinoma)
crc_all_stool <- read_csv("data/process/tables/alpha_RR_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_RR_ind_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V3", 3), rep("V4", 3), rep("V3-V4", 3), 
                    rep("V4", 6), rep("V3-V5", 3), rep("V3-V4", 3)))

crc_all_tissue <- read_csv("data/process/tables/alpha_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/alpha_RR_ind_tissue_results.csv")) %>% 
  mutate(region = c(rep("combined", 3), rep("V3-V5", 3), rep("V1-V2", 6), rep("V5-V6", 3), 
                    rep("V3-V4", 3), rep("V1-V3", 3)))


#### zeller = V4 - #CDB38B
#### hale = V3-V5 - #FFA500
#### brim = V1-V3 - #8B5A00
#### baxter = V4 - #CDB38B
#### weir = V4 - #CDB38B
#### wang = V3 - #B8860B
### flemer = V3-V4 - #E3CF57
### ahn = V3-V4 - #E3CF57
### lu = V3-V4 - #E3CF57
### dejea = V3-V5 - #FFA500
### geng = V1-V2 - #EEEE00
### sana = V1-V2 - #EEEE00
### burns = V5-V6 - #C76114
### chen = V1-V3 - #8B5A00

##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

adn_stool_graph <- adn_all_stool %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "zeller", "hale", "brim", "baxter"), 
                        labels = c( "Combined", "Zeller", "Hale", "Brim", "Baxter"))) %>% 
  filter(measure == "shannon") %>% 
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-2.2, 2.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("A") + 
  scale_color_manual(values = c( '#000000', '#8B5A00', '#FFA500', '#CDB38B')) + 
  annotate("text", label = paste("Adenoma (Stool)"), x = -1.8, y = 5.5, size = 1.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
  


crc_stool_graph <- crc_all_stool %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "zeller", "weir", "wang", "hale", 
                                   "flemer", "baxter", "ahn"), 
                        labels = c( "Combined", "Zeller", "Weir", "Wang",  "Hale", 
                                    "Flemer", "Baxter", "Ahn")), 
         region = factor(region, 
                         levels = c("combined", "V3", "V4", "V3-V4", "V3-V5"))) %>%  
  filter(measure == "shannon") %>% 
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-2.2, 2.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("B") + 
  scale_color_manual(values = c( '#000000', '#B8860B', '#CDB38B', '#E3CF57', '#FFA500')) + 
  annotate("text", label = paste("Carcinoma (Stool)"), x = -1.8, y = 8.4, size = 1.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


adn_tissue_graph <- adn_all_tissue %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "flemer", "lu"), 
                        labels = c( "Combined", "Flemer", "Lu"))) %>% 
  filter(measure == "shannon") %>% 
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-2.2, 3.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("C") + 
  scale_color_manual(values = c('#000000', '#E3CF57')) + 
  annotate("text", label = paste("Adenoma (Tissue)"), x = -1.65, y = 3.55, size = 1.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


crc_tissue_graph <- crc_all_tissue %>% 
  mutate(study = factor(study, 
                        levels = c("composite", "sana", "geng", "flemer", "dejea", "chen", "burns"), 
                        labels = c( "Combined", "Sanapareddy", "Geng", "Flemer", "Dejea", "Chen", "Burns")), 
         region = factor(region, 
                         levels = c("combined", "V1-V2", "V1-V3", "V3-V4", "V3-V5", "V5-V6"))) %>%  
  filter(measure == "shannon") %>% 
  ggplot(aes(log2(est), study, xmax=log2(upper), xmin=log2(lower), colour=region)) + 
  coord_cartesian(xlim=c(-4.2, 4.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("D") + 
  scale_color_manual(values = c( '#000000', '#EEEE00', '#8B5A00', '#E3CF57', '#FFA500', '#C76114')) + 
  annotate("text", label = paste("Carcinoma (Tissue)"), x = -3.25, y = 7.45, size = 1.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

shannon_RR <- grid.arrange(adn_stool_graph, crc_stool_graph, adn_tissue_graph, crc_tissue_graph)

ggsave("results/figures/shannon_RR.pdf", shannon_RR, width = 6, height = 6, dpi = 300)










