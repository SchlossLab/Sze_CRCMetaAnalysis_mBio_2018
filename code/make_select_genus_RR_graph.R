### Code to create Specific Genus RR pooled graphs
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse", "gridExtra", "viridis"))

# Load needed data tables (adenoma)
adn_all_stool <- read_csv("data/process/tables/adn_select_genus_RR_stool_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% arrange(pvalue, est) %>% 
  filter(est > 1) %>% slice(1:5) %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_RR_stool_composite.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
              mutate(study = "composite") %>% arrange(pvalue, est) %>% 
              filter(est < 1) %>% slice(1:5))

adn_ind_stool <- read_csv("data/process/tables/adn_select_genus_RR_stool_ind_results.csv") %>% 
  filter(measure %in% as.data.frame(adn_all_stool)[, "measure"])

adn_all_tissue <- read_csv("data/process/tables/adn_select_genus_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  arrange(pvalue, est) %>% 
  filter(est > 1) %>% slice(1:5) %>% 
  bind_rows(read_csv("data/process/tables/adn_select_genus_RR_tissue_composite.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
              mutate(study = "composite") %>% arrange(pvalue, est) %>% 
              filter(est < 1) %>% slice(1:5))

adn_ind_all_tissue <- read_csv("data/process/tables/adn_select_genus_RR_tissue_ind_results.csv") %>% 
  filter(measure %in% as.data.frame(adn_all_tissue)[, "measure"])


# Load in needed data tables (carcinoma)
crc_all_stool <- read_csv("data/process/tables/select_genus_RR_stool_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  arrange(pvalue, est) %>% 
  filter(est > 1) %>% slice(1:5) %>% 
  bind_rows(read_csv("data/process/tables/select_genus_RR_stool_composite.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
              mutate(study = "composite") %>% arrange(pvalue, est) %>% 
              filter(est < 1) %>% slice(1:5))

crc_ind_stool <- read_csv("data/process/tables/select_genus_RR_stool_ind_results.csv") %>% 
  filter(measure %in% as.data.frame(crc_all_stool)[, "measure"])


crc_all_tissue <- read_csv("data/process/tables/select_genus_RR_tissue_composite.csv") %>% 
  rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
  mutate(study = "composite") %>% 
  arrange(pvalue, est) %>% 
  filter(est > 1) %>% slice(1:5) %>% 
  bind_rows(read_csv("data/process/tables/select_genus_RR_tissue_composite.csv") %>% 
              rename(est = rr, lower = ci_lb, upper = ci_ub) %>% 
              mutate(study = "composite") %>% arrange(pvalue, est) %>% 
              filter(est < 1) %>% slice(1:5))

crc_ind_all_tissue <- read_csv("data/process/tables/select_genus_RR_tissue_ind_results.csv") %>% 
  filter(measure %in% as.data.frame(crc_all_tissue)[, "measure"])



##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

adn_stool_graph <- adn_all_stool %>% 
  mutate(measure = factor(measure, 
                        levels = c("Pyramidobacter", "Clostridium_XlVb", 
                                   "Candidatus_Saccharibacteria_unclassified", "Novosphingobium", 
                                   "Bacteroidales_unclassified", "Lachnospiraceae_unclassified", 
                                   "Lactococcus", "Clostridium_XI", "Firmicutes_unclassified", 
                                   "Clostridiaceae_1_unclassified"), 
                        labels = c("Pyramidobacter", "Clostridium XlVb", 
                                   "Candidatus\nSaccharibacteria", "Novosphingobium", 
                                   "Bacteroidales", "Lachnospiraceae", 
                                   "Lactococcus", "Clostridium XI", "Firmicutes", 
                                   "Clostridiaceae 1"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-2.5, 2.5)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 2.5, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("A") + 
  scale_color_manual(values = c('#B0171F', '#B0171F', '#B0171F', '#B0171F', '#B0171F', 
                                '#0000EE', '#0000EE', '#0000EE', '#0000EE', '#0000EE')) + 
  scale_y_discrete(labels=expression(
    italic(Pyramidobacter), italic(Clostridium~XlVb), 
    italic(Candidatus~Saccharibacteria), 
    italic(Novosphingobium), italic(Bacteroidales), italic(Lachnospiraceae), 
    italic(Lactococcus), italic(Clostridium~XI), italic(Firmicutes), 
    italic(Clostridiaceae)~1)) +
  annotate("text", label = paste("Adenoma\n(Stool)"), x = 2, y = 9.8, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.9, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8))


crc_stool_graph <- crc_all_stool %>% 
  mutate(measure = factor(measure, 
                          levels = c("Peptostreptococcus", "Porphyromonas", "Parvimonas", 
                                     "Fusobacterium", "Escherichia.Shigella", "Roseburia", 
                                     "Ruminococcus", "Lachnospiraceae_unclassified", 
                                     "Clostridium_XI", "Clostridiaceae_1_unclassified"), 
                          labels = c("Peptostreptococcus", "Porphyromonas", "Parvimonas", 
                                     "Fusobacterium", "Escherichia/Shigella", "Roseburia", 
                                     "Ruminococcus", "Lachnospiraceae", 
                                     "Clostridium XI", "Clostridiaceae 1"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-2.5, 2.5)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 2.5, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("B") + 
  scale_color_manual(values = c('#B0171F', '#B0171F', '#B0171F', '#B0171F', '#B0171F', 
                                '#0000EE', '#0000EE', '#0000EE', '#0000EE', '#0000EE')) + 
  scale_y_discrete(labels=expression(
    italic(Peptostreptococcus), italic(Porphyromonas), italic(Parvimonas), italic(Fusobacterium), 
    italic(Escherichia)/italic(Shigella), italic(Roseburia), italic(Ruminococcus), 
    italic(Lachnospiraceae), italic(Clostridium~XI), italic(Clostridiaceae)~1)) +
  annotate("text", label = paste("Carcinoma\n(Stool)"), x = 2, y = 9.8, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.6, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8))


adn_tissue_graph <- adn_all_tissue %>% 
  mutate(measure = factor(measure, 
                          levels = c("Selenomonas", "Enterobacter", "Rothia", 
                                     "Micrococcaceae_unclassified", "Achromobacter", 
                                     "Lachnospiraceae_unclassified", "Butyricicoccus", 
                                     "Clostridiaceae_1_unclassified", "Parasutterella", 
                                     "Pseudoflavonifractor"), 
                          labels = c("Selenomonas", "Enterobacter", "Rothia", 
                                     "Micrococcaceae", "Achromobacter", 
                                     "Lachnospiraceae", "Butyricicoccus", 
                                     "Clostridiaceae 1", "Parasutterella", 
                                     "Pseudoflavonifractor"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-2.5, 2.5)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 2.5, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("C") + 
  scale_color_manual(values = c('#B0171F', '#B0171F', '#B0171F', '#B0171F', '#B0171F', 
                                '#0000EE', '#0000EE', '#0000EE', '#0000EE', '#0000EE')) + 
  scale_y_discrete(labels=expression(
    italic(Selenomonas), italic(Enterobacter), italic(Rothia), italic(Micrococcaceae), 
    italic(Achromobacter), italic(Lachnospiraceae), italic(Butyricicoccus), 
    italic(Clostridiaceae)~1, italic(Parasutterella), italic(Pseudoflavonifractor))) + 
  annotate("text", label = paste("Adenoma\n(Tissue)"), x = 2, y = 9.8, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.6, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8))


crc_tissue_graph <- crc_all_tissue %>% 
  mutate(measure = factor(measure, 
                          levels = c("Campylobacter", "Leptotrichia", "Lactobacillus", 
                                     "Anaerococcus", "Fusobacterium", "Bacteroides", 
                                     "Corynebacterium", "Blautia", "Dorea", "Clostridium_XlVb"), 
                          labels = c("Campylobacter", "Leptotrichia", "Lactobacillus", 
                                     "Anaerococcus", "Fusobacterium", "Bacteroides", 
                                     "Corynebacterium", "Blautia", "Dorea", "Clostridium XlVb"))) %>% 
  ggplot(aes(log2(est), measure, xmax=log2(upper), xmin=log2(lower), colour=measure)) + 
  coord_cartesian(xlim=c(-2.5, 2.5)) + 
  geom_vline(xintercept = 0.0, linetype=2, alpha=0.75) + 
  geom_errorbarh(alpha=0.5, size = 1, height=0, show.legend = F) + 
  geom_point(size = 2.5, show.legend = F) + 
  labs(x = expression(Log["2"]~Relative~Risk), y = "") + theme_bw() + ggtitle("D") + 
  scale_color_manual(values = c('#B0171F', '#B0171F', '#B0171F', '#B0171F', '#B0171F', 
                                '#0000EE', '#0000EE', '#0000EE', '#0000EE', '#0000EE')) + 
  scale_y_discrete(labels=expression(
    italic(Campylobacter), italic(Leptotrichia), italic(Lactobacillus), italic(Anaerococcus), 
    italic(Fusobacterium), italic(Bacteroides), italic(Corynebacterium), italic(Blautia), 
    italic(Dorea), italic(Clostridium~XlVb))) + 
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

