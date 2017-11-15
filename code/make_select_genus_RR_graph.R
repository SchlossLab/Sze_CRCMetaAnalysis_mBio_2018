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
  bind_rows(read_csv("data/process/tables/adn_select_genus_RR_stool_ind_results.csv")) %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_inc_4_stool.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub)) %>% 
  mutate(study = ifelse(is.na(study), invisible("composite"), invisible(study)))

adn_all_tissue <- read_csv("data/process/tables/adn_select_genus_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_RR_tissue_ind_results.csv")) %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_inc_4_tissue.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub)) %>% 
  mutate(study = ifelse(is.na(study), invisible("composite"), invisible(study)))

# Load in needed data tables (carcinoma)
crc_all_stool <- read_csv("data/process/tables/select_genus_RR_stool_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/select_genus_RR_stool_ind_results.csv")) %>% 
  bind_rows(read_csv("data/process/tables/select_genus_inc_4_stool.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub)) %>% 
  mutate(study = ifelse(is.na(study), invisible("composite"), invisible(study))) %>% 
  filter(!(measure %in% c("one_or_more", "two_or_more", "three_or_more", "four_only")))


crc_all_tissue <- read_csv("data/process/tables/select_genus_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  bind_rows(read_csv("data/process/tables/select_genus_RR_tissue_ind_results.csv")) %>% 
  bind_rows(read_csv("data/process/tables/select_genus_inc_4_tissue.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub)) %>% 
  mutate(study = ifelse(is.na(study), invisible("composite"), invisible(study))) %>% 
  filter(!(measure %in% c("one_or_more", "two_or_more", "three_or_more", "four_only")))


##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

adn_stool_graph <- adn_all_stool %>% 
  filter(study == "composite", 
         measure != "all_four", 
         measure != "total_four") %>% 
  mutate(measure = factor(measure, 
                        levels = c("four_counts", "three_counts", "two_counts", "one_counts", 
                                   "Fusobacterium", "Porphyromonas", "Peptostreptococcus", "Parvimonas"), 
                        labels = c("All 4", "At Least 3", "At Least 2", "At Least 1", "Fusobacterium", 
                                   "Porphyromonas", "Peptostreptococcus", "Parvimonas"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-4.2, 4.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("A") + 
  scale_color_manual(values = c('#000000', '#000000', '#000000', '#000000', 
                                '#778899', '#778899', '#778899', '#778899')) + 
  scale_y_discrete(labels=expression(All~4, At~Least~3, At~Least~2, At~Least~1, italic(Fusobacterium), 
                                     italic(Porphyromonas), italic(Peptostreptococcus), italic(Parvimonas))) +
  annotate("text", label = paste("Adenoma (Stool)"), x = -2.9, y = 8.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.8, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_stool_graph <- crc_all_stool %>% 
  filter(study == "composite", 
         measure != "all_four", 
         measure != "total_four") %>% 
  mutate(measure = factor(measure, 
                          levels = c("four_counts", "three_counts", "two_counts", "one_counts", 
                                     "Fusobacterium", "Porphyromonas", "Peptostreptococcus", "Parvimonas"), 
                          labels = c("All 4", "At Least 3", "At Least 2", "At Least 1", "Fusobacterium", 
                                     "Porphyromonas", "Peptostreptococcus", "Parvimonas"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-4.2, 4.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("B") + 
  scale_color_manual(values = c('#000000', '#000000', '#000000', '#000000', 
                                '#778899', '#778899', '#778899', '#778899')) + 
  scale_y_discrete(labels=expression(All~4, At~Least~3, At~Least~2, At~Least~1, italic(Fusobacterium), 
                                     italic(Porphyromonas), italic(Peptostreptococcus), italic(Parvimonas))) +
  annotate("text", label = paste("Carcinoma (Stool)"), x = -2.8, y = 8.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.8, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


adn_tissue_graph <- adn_all_tissue %>% 
  filter(study == "composite", 
         measure != "all_four", 
         measure != "total_four") %>% 
  mutate(measure = factor(measure, 
                          levels = c("four_counts", "three_counts", "two_counts", "one_counts", 
                                     "Fusobacterium", "Porphyromonas", "Peptostreptococcus", "Parvimonas"), 
                          labels = c("All 4", "At Least 3", "At Least 2", "At Least 1", "Fusobacterium", 
                                     "Porphyromonas", "Peptostreptococcus", "Parvimonas"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-4.2, 4.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("C") + 
  scale_color_manual(values = c('#000000', '#000000', '#000000', '#000000', 
                                '#778899', '#778899', '#778899', '#778899')) + 
  scale_y_discrete(labels=expression(All~4, At~Least~3, At~Least~2, At~Least~1, italic(Fusobacterium), 
                                     italic(Porphyromonas), italic(Peptostreptococcus), italic(Parvimonas))) + 
  annotate("text", label = paste("Adenoma (Tissue)"), x = -2.8, y = 8.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.8, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


crc_tissue_graph <- crc_all_tissue %>% 
  filter(study == "composite", 
         measure != "all_four", 
         measure != "total_four") %>% 
  mutate(measure = factor(measure, 
                          levels = c("four_counts", "three_counts", "two_counts", "one_counts", 
                                     "Fusobacterium", "Porphyromonas", "Peptostreptococcus", "Parvimonas"), 
                          labels = c("All 4", "At Least 3", "At Least 2", "At Least 1", "Fusobacterium", 
                                     "Porphyromonas", "Peptostreptococcus", "Parvimonas"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-4.2, 4.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 3, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("D") + 
  scale_color_manual(values = c('#000000', '#000000', '#000000', '#000000', 
                                '#778899', '#778899', '#778899', '#778899')) + 
  scale_y_discrete(labels=expression(All~4, At~Least~3, At~Least~2, At~Least~1, italic(Fusobacterium), 
                                     italic(Porphyromonas), italic(Peptostreptococcus), italic(Parvimonas))) + 
  annotate("text", label = paste("Carcinoma (Tissue)"), x = -2.7, y = 8.4, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.8, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

select_genus_combined <- grid.arrange(adn_stool_graph, crc_stool_graph, 
                                      adn_tissue_graph, crc_tissue_graph)

ggsave("results/figures/Figure3.pdf", 
       select_genus_combined, width = 8, height = 6, dpi = 300)

