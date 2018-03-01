### Code for the crc most imp taxa in select models
### Marc Sze

# Load in needed generic functions
source("code/functions.R")

# Load in needed libraries
loadLibs(c("tidyverse"))


ind_crc_stool_auc <- read_csv("data/process/tables/ind_genera_auc_stool.csv") %>%  
  mutate(auc = ifelse(taxa == "Ruminococcus" | taxa == "Clostridium_XI", invisible(1-auc), invisible(auc)))

ind_unmatched_tissue_crc_auc <- read_csv("data/process/tables/ind_genera_auc_unmatched_tissue.csv") %>% 
  mutate(auc = ifelse(taxa == "Dorea" | taxa == "Blautia", invisible(1-auc), invisible(auc)))


crc_tissue_model_aucs <- read_csv("data/process/tables/ALL_genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv") %>% 
  mutate(type = "unmatched") %>% 
  bind_rows(read_csv("data/process/tables/ALL_genus_matched_tissue_RF_fullvsselect_pvalue_summary.csv") %>% 
              mutate(type = "matched")) %>% 
  gather(key = model_type, value = AUC, full_model, select_model) %>% 
  mutate(AUC = as.numeric(AUC)) %>% 
  bind_rows(read_csv("data/process/tables/matched_tissue_rf_otu_random_comparison_summary.csv") %>% 
              mutate(type = "matched", 
                     model_type = "full_otu", 
                     train_model = study, 
                     act_mean_auc = as.numeric(act_mean_auc)) %>% 
              rename(AUC = act_mean_auc) %>% 
              select(pvalue, study, train_model, type, model_type, AUC) %>% 
              bind_rows(read_csv("data/process/tables/unmatched_tissue_rf_otu_random_comparison_summary.csv") %>% 
                          mutate(type = "unmatched", 
                                 model_type = "full_otu", 
                                 train_model = study, 
                                 act_mean_auc = as.numeric(act_mean_auc)) %>% 
                          rename(AUC = act_mean_auc) %>% 
                          select(pvalue, study, train_model, type, model_type, AUC)))



crc_stool_model_aucs <- read_csv("data/process/tables/ALL_genus_stool_RF_fullvsselect_pvalue_summary.csv") %>% 
  gather(key = model_type, value = AUC, full_model, select_model) %>% 
  mutate(AUC = as.numeric(AUC), pvalue = as.numeric(pvalue)) %>% 
  bind_rows(read_csv("data/process/tables/stool_rf_otu_random_comparison_summary.csv") %>% 
              mutate(model_type = "full_otu", 
                     train_model = study, 
                     act_mean_auc = as.numeric(act_mean_auc)) %>% 
              rename(AUC = act_mean_auc) %>% 
              select(pvalue, study, train_model, model_type, AUC))


##############################################################################################
############################## List of code to make figures ##################################
##############################################################################################

# Function to find pvalue
get_pvalue <- function(ind_taxa, modelType, full_taxa_model_data, ind_taxa_model_data){
  
  
  temp_pvalue <- try(t.test(as.data.frame(filter(full_taxa_model_data, 
                                             model_type == "full_model", study == train_model))[, "AUC"], 
                        as.data.frame(filter(ind_taxa_model_data, taxa == ind_taxa))[, "auc"], 
                        paired = TRUE)$p.value)
  
  
  return(temp_pvalue)
  
}



# Function to find tissue pvalue
get_tissue_pvalue <- function(ind_taxa, modelType, full_taxa_model_data, ind_taxa_model_data){
  
  
  temp_pvalue <- try(t.test(as.data.frame(filter(full_taxa_model_data, 
                                                 model_type == modelType, type == "unmatched", 
                                                 study == train_model))[, "AUC"], 
                            as.data.frame(filter(ind_taxa_model_data, taxa == ind_taxa))[, "auc"], 
                            paired = TRUE)$p.value)
  
  
  return(temp_pvalue)
  
}



##############################################################################################
############### Run the actual programs to make the figure ###################################
##############################################################################################

# get taxa of interest for stool
taxa_of_int <- unique(ind_crc_stool_auc$taxa)

# get taxa of interest for unmatched tissue
taxa_of_int_unmatched <- unique(ind_unmatched_tissue_crc_auc$taxa)

# get pvalues comparisons of stool
stool_tests <- sapply(taxa_of_int, 
                      function(x) get_pvalue(x, "full_model", crc_stool_model_aucs, ind_crc_stool_auc), simplify = T) %>% 
  bind_cols(taxa = taxa_of_int, pvalue = .) %>% 
  mutate(bh = p.adjust(pvalue, method = "BH"))


unmatched_tissue_tests <- sapply(taxa_of_int_unmatched, 
                      function(x) get_tissue_pvalue(x, "full_model", crc_tissue_model_aucs, 
                                             ind_unmatched_tissue_crc_auc), simplify = T) %>% 
  bind_cols(taxa = taxa_of_int_unmatched, pvalue = .) %>% 
  mutate(bh = p.adjust(pvalue, method = "BH"))

# Write out the data
write_csv(stool_tests, "data/process/tables/stool_ind_vs_full_taxa_pvalue_summary.csv")

write_csv(unmatched_tissue_tests, "data/process/tables/unmatched_tissue_ind_vs_full_taxa_pvalue_summary.csv")






