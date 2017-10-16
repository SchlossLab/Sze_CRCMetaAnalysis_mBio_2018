### Code to create Specific Genus RR pooled graphs
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))

# Load needed data tables (adenoma)
adn_all_stool <- read_csv("data/process/tables/adn_select_genus_RR_stool_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_RR_stool_ind_results.csv"))

adn_all_tissue <- read_csv("data/process/tables/adn_select_genus_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_RR_tissue_ind_results.csv"))

# Load in needed data tables (carcinoma)
crc_all_stool <- read_csv("data/process/tables/select_genus_RR_stool_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/select_genus_RR_stool_ind_results.csv"))


crc_all_tissue <- read_csv("data/process/tables/select_genus_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/select_genus_RR_tissue_ind_results.csv"))


##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

adn_stool_graph <- adn_all_stool %>% 
  filter(study == "composite") %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-4.2, 4.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("A") + 
  scale_color_manual(values = c('#8B5A00', '#FFA500', '#FFC125', '#E3CF57')) + 
  annotate("text", label = paste("Adenoma (Stool)"), x = -2.8, y = 4.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.6, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_stool_graph <- crc_all_stool %>% 
  filter(study == "composite") %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-4.2, 4.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("B") + 
  scale_color_manual(values = c('#B0171F', '#FF4040', '#CD5C5C', '#FFC1C1')) + 
  annotate("text", label = paste("Carcinoma (Stool)"), x = -2.8, y = 4.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.6, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


adn_tissue_graph <- adn_all_tissue %>% 
  filter(study == "composite") %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-4.2, 4.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("C") + 
  scale_color_manual(values = c('#8B5A00', '#FFA500', '#FFC125', '#E3CF57')) + 
  annotate("text", label = paste("Adenoma (Tissue)"), x = -2.8, y = 4.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.6, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_tissue_graph <- crc_all_tissue %>% 
  filter(study == "composite") %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-4.2, 4.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("D") + 
  scale_color_manual(values = c('#B0171F', '#FF4040', '#CD5C5C', '#FFC1C1')) + 
  annotate("text", label = paste("Carcinoma (Tissue)"), x = -2.8, y = 4.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.6, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

select_genus_combined <- grid.arrange(adn_stool_graph, crc_stool_graph, 
                                      adn_tissue_graph, crc_tissue_graph)

ggsave("results/figures/select_genus_RR_composite.pdf", 
       select_genus_combined, width = 8, height = 6, dpi = 300)

