### Code to create Specific Genus RR pooled graphs
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))

# Load needed data tables (adenoma)
adn_all_stool <- read_csv("data/process/tables/adn_select_genus_OR_stool_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% arrange(pvalue, est) %>% 
  filter(est > 1) %>% slice(1:5) %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_OR_stool_composite.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
              mutate(study = "composite") %>% arrange(pvalue, est) %>% 
              filter(est < 1) %>% slice(1:5)) %>% 
  mutate(high_low = ifelse(est > 1, invisible("high"), invisible("low")))

adn_ind_stool <- read_csv("data/process/tables/adn_select_genus_OR_stool_ind_results.csv") %>% 
  filter(measure %in% as.data.frame(adn_all_stool)[, "measure"]) %>% 
  bind_rows(adn_all_stool) %>% 
  mutate(high_low = ifelse(measure %in% c("Clostridium_XlVb", "Porphyromonas", 
                                          "Novosphingobium", "Bacteroidales_unclassified", 
                                          "Catenibacterium"), 
                           invisible("high"), invisible("low")))


adn_all_tissue <- read_csv("data/process/tables/adn_select_genus_OR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  arrange(pvalue, est) %>% 
  filter(est > 1) %>% slice(1:5) %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_OR_tissue_composite.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
              mutate(study = "composite") %>% arrange(pvalue, est) %>% 
              filter(est < 1) %>% slice(1:5)) %>% 
  mutate(high_low = ifelse(est > 1, invisible("high"), invisible("low")))

adn_ind_all_tissue <- read_csv("data/process/tables/adn_select_genus_OR_tissue_ind_results.csv") %>% 
  filter(measure %in% as.data.frame(adn_all_tissue)[, "measure"]) %>% 
  mutate(high_low = ifelse(measure %in% c("Pseudomonas", "Howardella", 
                                          "Rothia", "Enterobacter", 
                                          "Puniceicoccaceae_unclassified"), 
                           invisible("high"), invisible("low")))


# Load in needed data tables (carcinoma)
crc_all_stool <- read_csv("data/process/tables/select_genus_OR_stool_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  arrange(pvalue, est) %>% 
  filter(est > 1) %>% slice(1:5) %>% 
  bind_rows(read_csv("data/process/tables/select_genus_OR_stool_composite.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
              mutate(study = "composite") %>% arrange(pvalue, est) %>% 
              filter(est < 1) %>% slice(1:5)) %>% 
  mutate(high_low = ifelse(est > 1, invisible("high"), invisible("low")))

crc_ind_stool <- read_csv("data/process/tables/select_genus_OR_stool_ind_results.csv") %>% 
  filter(measure %in% as.data.frame(crc_all_stool)[, "measure"]) %>% 
  mutate(high_low = ifelse(measure %in% c("Porphyromonas", "Peptostreptococcus", 
                                          "Parvimonas", 
                                          "Fusobacterium", "Escherichia.Shigella"), 
                           invisible("high"), invisible("low")))


crc_all_tissue <- read_csv("data/process/tables/select_genus_OR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  arrange(pvalue, est) %>% 
  filter(est > 1) %>% slice(1:5) %>% 
  bind_rows(read_csv("data/process/tables/select_genus_OR_tissue_composite.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
              mutate(study = "composite") %>% arrange(pvalue, est) %>% 
              filter(est < 1) %>% slice(1:5)) %>% 
  mutate(high_low = ifelse(est > 1, invisible("high"), invisible("low")))

crc_ind_all_tissue <- read_csv("data/process/tables/select_genus_OR_tissue_ind_results.csv") %>% 
  filter(measure %in% as.data.frame(crc_all_tissue)[, "measure"]) %>% 
  mutate(high_low = ifelse(measure %in% c("Clostridium_sensu_stricto", "Campylobacter", 
                                          "Leptotrichia", "Fusobacterium", 
                                          "Parvimonas"), 
                           invisible("high"), invisible("low")))



##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

adn_stool_graph <- adn_all_stool %>% 
  mutate(measure = factor(measure, 
                        levels = c("Clostridium_XlVb", "Porphyromonas", 
                                   "Novosphingobium", "Bacteroidales_unclassified", 
                                   "Catenibacterium", "Lachnospiraceae_unclassified", 
                                   "Clostridium_XI", "Clostridiaceae_1_unclassified", 
                                   "Lactococcus", "Lactobacillales_unclassified"), 
                        labels = c("Clostridium XlVb", "Porphyromonas", 
                                   "Novosphingobium", "Bacteroidales", 
                                   "Catenibacterium",  "Lachnospiraceae", 
                                   "Clostridium XI", "Clostridiaceae 1", 
                                   "Lactococcus", "Lactobacillales"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), 
             xmin=log2(lower), color = high_low)) + 
  coord_cartesian(xlim=c(-5.2, 5.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(aes(group = high_low), size = 2.5, show.legend = F) + 
  geom_point(data = adn_ind_stool, aes(log2(est), 
    factor(measure, 
           levels = c("Clostridium_XlVb", "Porphyromonas", 
                      "Novosphingobium", "Bacteroidales_unclassified", 
                      "Catenibacterium", "Lachnospiraceae_unclassified", 
                      "Clostridium_XI", "Clostridiaceae_1_unclassified", 
                      "Lactococcus", "Lactobacillales_unclassified"), 
           labels = c("Clostridium XlVb", "Porphyromonas", 
                      "Novosphingobium", "Bacteroidales", 
                      "Catenibacterium",  "Lachnospiraceae", 
                      "Clostridium XI", "Clostridiaceae 1", 
                      "Lactococcus", "Lactobacillales")), group = high_low, color = high_low), 
             show.legend = F, alpha = 0.5, size = 1.25) + 
  labs(x = expression(Log["2"]~Odds~Ratio), y = "") + theme_bw() + ggtitle("A") + 
  scale_color_manual(values = c('#B0171F', '#0000EE')) + 
  scale_y_discrete(labels=expression(
    italic(Clostridium~XlVb), italic(Porphyromonas), 
    italic(Novosphingobium), 
    italic(Bacteroidales), italic(Catenibacterium), italic(Lachnospiraceae), 
    italic(Clostridium~XI), italic(Clostridiaceae)~1, italic(Lactococcus), 
    italic(Lactobacillales))) +
  annotate("text", label = paste("Adenoma\n(Stool)"), x = 4.5, y = 9.8, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.9, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8))


crc_stool_graph <- crc_all_stool %>% 
  mutate(measure = factor(measure, 
                          levels = c("Porphyromonas", "Peptostreptococcus", "Parvimonas", 
                                     "Fusobacterium", "Escherichia.Shigella", "Ruminococcus", 
                                     "Clostridium_XI", "Roseburia", "Clostridiaceae_1_unclassified"
                                     , "Lachnospiraceae_unclassified"), 
                          labels = c("Porphyromonas", "Peptostreptococcus", "Parvimonas", 
                                     "Fusobacterium", "Escherichia/Shigella", "Ruminococcus", 
                                     "Clostridium XI", "Roseburia", "Clostridiaceae 1", 
                                     "Lachnospiraceae"))) %>% 
  ggplot(aes(log2(est), measure, 
             xmax=log2(upper), xmin=log2(lower), colour=high_low)) + 
  coord_cartesian(xlim=c(-5.2, 5.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(aes(group = high_low), size = 2.5, show.legend = F) + 
  geom_point(data = crc_ind_stool, aes(log2(est), 
        factor(measure, 
               levels = c("Porphyromonas", "Peptostreptococcus", "Parvimonas", 
                          "Fusobacterium", "Escherichia.Shigella", "Ruminococcus", 
                          "Clostridium_XI", "Roseburia", "Clostridiaceae_1_unclassified", 
                          "Lachnospiraceae_unclassified"), 
               labels = c("Porphyromonas", "Peptostreptococcus", "Parvimonas", 
                          "Fusobacterium", "Escherichia/Shigella", "Ruminococcus", 
                          "Clostridium XI", "Roseburia", "Clostridiaceae 1", 
                          "Lachnospiraceae")), 
        group = high_low, color = high_low), 
             show.legend = F, alpha = 0.5, size = 1.25) + 
  labs(x = expression(Log["2"]~Odds~Ratio), y = "") + theme_bw() + ggtitle("B") + 
  scale_color_manual(values = c('#B0171F', '#0000EE')) + 
  scale_y_discrete(labels=expression(
    italic(Porphyromonas), italic(Peptostreptococcus), italic(Parvimonas), italic(Fusobacterium), 
    italic(Escherichia)/italic(Shigella), italic(Ruminococcus), italic(Clostridium~XI), 
    italic(Roseburia), italic(Clostridiaceae)~1, italic(Lachnospiraceae))) +
  annotate("text", label = paste("Carcinoma\n(Stool)"), x = 4.5, y = 9.8, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.6, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8))


adn_tissue_graph <- adn_all_tissue %>% 
  mutate(measure = factor(measure, 
                          levels = c("Pseudomonas", "Howardella", "Rothia", 
                                     "Enterobacter", "Puniceicoccaceae_unclassified", 
                                     "Lachnospiraceae_unclassified", "Blautia", 
                                     "Anaerostipes", "Butyricicoccus", 
                                     "Parasutterella"), 
                          labels = c("Pseudomonas", "Howardella", "Rothia", 
                                     "Enterobacter", "Puniceicoccaceae", 
                                     "Lachnospiraceae", "Blautia", 
                                     "Anaerostipes", "Butyricicoccus", 
                                     "Parasutterella"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=high_low)) + 
  coord_cartesian(xlim=c(-5.2, 5.2)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(aes(group = high_low), size = 2.5, show.legend = F) + 
  geom_point(data = adn_ind_all_tissue, aes(log2(est), 
              factor(measure, 
                     levels = c("Pseudomonas", "Howardella", "Rothia", 
                                "Enterobacter", "Puniceicoccaceae_unclassified", 
                                "Lachnospiraceae_unclassified", "Blautia", 
                                "Anaerostipes", "Butyricicoccus", 
                                "Parasutterella"), 
                     labels = c("Pseudomonas", "Howardella", "Rothia", 
                                "Enterobacter", "Puniceicoccaceae", 
                                "Lachnospiraceae", "Blautia", 
                                "Anaerostipes", "Butyricicoccus", 
                                "Parasutterella")), 
                                       group = high_low, color = high_low), 
             show.legend = F, alpha = 0.5, size = 1.25) + 
  labs(x = expression(Log["2"]~Odds~Ratio), y = "") + theme_bw() + ggtitle("C") + 
  scale_color_manual(values = c('#B0171F', '#0000EE')) + 
  scale_y_discrete(labels=expression(
    italic(Pseudomonas), italic(Howardella), italic(Rothia), italic(Enterobacter), 
    italic(Puniceicoccaceae), italic(Lachnospiraceae), italic(Blautia), 
    italic(Anaerostipes), italic(Butyricicoccus), italic(Parasutterella))) + 
  annotate("text", label = paste("Adenoma\n(Tissue)"), x = 4.5, y = 9.8, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.6, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8))


crc_tissue_graph <- crc_all_tissue %>% 
  mutate(measure = factor(measure, 
                          levels = c("Clostridium_sensu_stricto", "Campylobacter", "Leptotrichia", 
                                     "Fusobacterium", "Parvimonas", "Blautia", "Corynebacterium", 
                                     "Bacteroides", "Clostridium_XlVb", "Ruminococcus2"), 
                          labels = c("Clostridium sensu stricto", "Campylobacter", "Leptotrichia", 
                                     "Fusobacterium", "Parvimonas", "Blautia", "Corynebacterium", 
                                     "Bacteroides", "Clostridium XlVb", "Ruminococcus"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=high_low)) + 
  coord_cartesian(xlim=c(-2.5, 2.5)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(aes(group = high_low), size = 2.5, show.legend = F) + 
  geom_point(data = crc_ind_all_tissue, aes(log2(est), 
                factor(measure, 
                       levels = c("Clostridium_sensu_stricto", "Campylobacter", "Leptotrichia", 
                                  "Fusobacterium", "Parvimonas", "Blautia", "Corynebacterium", 
                                  "Bacteroides", "Clostridium_XlVb", "Ruminococcus2"), 
                       labels = c("Clostridium sensu stricto", "Campylobacter", "Leptotrichia", 
                                  "Fusobacterium", "Parvimonas", "Blautia", "Corynebacterium", 
                                  "Bacteroides", "Clostridium XlVb", "Ruminococcus")), 
                                            group = high_low, color = high_low), 
             show.legend = F, alpha = 0.5, size = 1.25) +  
  labs(x = expression(Log["2"]~Odds~Ratio), y = "") + theme_bw() + ggtitle("D") + 
  scale_color_manual(values = c('#B0171F', '#0000EE')) + 
  scale_y_discrete(labels=expression(
    italic(Clostridium~sensu~stricto), italic(Campylobacter), italic(Leptotrichia), italic(Fusobacterium), 
    italic(Parvimonas), italic(Blautia), italic(Corynebacterium), italic(Bacteroides), 
    italic(Clostridium~XlVb), italic(Ruminococcus))) + 
  annotate("text", label = paste("Carcinoma\n(Tissue)"), x = 2, y = 9.8, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.5, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

select_genus_combined <- grid.arrange(adn_stool_graph, crc_stool_graph, 
                                      adn_tissue_graph, crc_tissue_graph)

ggsave("results/figures/Figure3.pdf", 
       select_genus_combined, width = 8, height = 6, dpi = 300)

