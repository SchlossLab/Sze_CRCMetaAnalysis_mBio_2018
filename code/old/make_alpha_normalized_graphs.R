### Combined Alpha graphs
### Graphs of z-score normalized and power transformed data
### Marc Sze


# Load in needed functions and libraries
source('code/functions.R')

# Load needed libraries
loadLibs(c("tidyverse", "gridExtra"))


# Load in needed data

stool_data <- read_csv("data/process/tables/stool_normalized_alpha_all.csv") %>% 
  select(group, disease, sobs, shannon, shannoneven) %>% 
  gather(key = alpha_metric, value = measurement, sobs, shannon, shannoneven)

matched_tissue_data <- read_csv("data/process/tables/matched_tissue_normalized_alpha_all_data.csv") %>% 
  select(id, disease, sobs, shannon, shannoneven) %>% 
  gather(key = alpha_metric, value = measurement, sobs, shannon, shannoneven)

unmatched_tissue_data <- read_csv("data/process/tables/unmatched_tissue_normalized_alpha_all_data.csv") %>% 
  select(id, disease, sobs, shannon, shannoneven) %>% 
  gather(key = alpha_metric, value = measurement, sobs, shannon, shannoneven)



##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

stool_graph <- stool_data %>% 
  mutate(disease = factor(disease, 
                          levels = c("control", "polyp", "cancer"), 
                          labels = c("Control", "Adenoma", "Carcinoma")), 
         alpha_metric = factor(alpha_metric, 
                               levels = c("sobs", "shannoneven", "shannon"), 
                               labels = c("Observed OTUs", "Evenness", "Shannon Diversity"))) %>% 
  ggplot(aes(disease, measurement, color = disease)) + 
  geom_boxplot(show.legend = F) + facet_grid(. ~ alpha_metric) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c('#228B22', '#FFD700', '#DC143C')) + 
  labs(x = "", y = "Z-Score") + theme_bw() + ggtitle("A") + 
  annotate("text", label = paste("Stool"), x = 2, y = 3.5, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))



unmatched_tissue_graph <- unmatched_tissue_data %>% 
  mutate(disease = factor(disease, 
                          levels = c("control", "polyp", "cancer"), 
                          labels = c("Control", "Adenoma", "Carcinoma")), 
         alpha_metric = factor(alpha_metric, 
                               levels = c("sobs", "shannoneven", "shannon"), 
                               labels = c("Observed OTUs", "Evenness", "Shannon Diversity"))) %>% 
  ggplot(aes(disease, measurement, color = disease)) + 
  geom_boxplot(show.legend = F) + facet_grid(. ~ alpha_metric) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c('#228B22', '#FFD700', '#DC143C')) + 
  labs(x = "", y = "Z-Score") + theme_bw() + ggtitle("B") + 
  annotate("text", label = paste("Unmatched Tissue"), x = 2, y = 3, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


matched_tissue_graph <- matched_tissue_data %>% 
  mutate(disease = factor(disease, 
                          levels = c("control", "polyp", "cancer"), 
                          labels = c("Control", "Adenoma", "Carcinoma")), 
         alpha_metric = factor(alpha_metric, 
                               levels = c("sobs", "shannoneven", "shannon"), 
                               labels = c("Observed OTUs", "Evenness", "Shannon Diversity"))) %>% 
  ggplot(aes(disease, measurement, color = disease)) + 
  geom_boxplot(show.legend = F) + facet_grid(. ~ alpha_metric) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+ 
  scale_color_manual(values = c('#228B22', '#FFD700', '#DC143C')) + 
  labs(x = "", y = "Z-Score") + theme_bw() + ggtitle("C") + 
  annotate("text", label = paste("Matched Tissue"), x = 2, y = 3, size = 2.5) + 
  theme(plot.title = element_text(face="bold", hjust = -0.07, size = 20), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 10))


##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################


alpha_all_graph <- grid.arrange(stool_graph, unmatched_tissue_graph, matched_tissue_graph)


ggsave("results/figures/Figure1.pdf", 
       alpha_all_graph, width = 6.5, height = 9, dpi = 300)



