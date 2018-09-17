STUDIES = sana brim geng weir dejea baxter lu ahn zeller burns wang chen flemer hale
REFS = data/references
FIGS = results/figures
TABLES = data/process/tables
PROC = data/process
FINAL = submission/


# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam ...
print-%:
	@echo '$*=$($*)'

################################################################################
#
# Part 1: Get the references
#
# We will need several reference files to complete the analyses including the
# SILVA reference alignment and RDP reference taxonomy.
#
################################################################################

# We want the latest greatest reference alignment and the SILVA reference
# alignment is the best reference alignment on the market. This version is from
# v123 and described at http://blog.mothur.org/2015/12/03/SILVA-v123-reference-files/
# We will use the SEED v. 123, which contain 12,083 bacterial sequences. This
# also contains the reference taxonomy. We will limit the databases to only
# include bacterial sequences.

$(REFS)/silva.seed.align :
	wget -N http://mothur.org/w/images/1/15/Silva.seed_v123.tgz
	tar xvzf Silva.seed_v123.tgz silva.seed_v123.align silva.seed_v123.tax
	mothur "#get.lineage(fasta=silva.seed_v123.align, taxonomy=silva.seed_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v123.pick.align, processors=8)"
	mv silva.seed_v123.pick.align $(REFS)/silva.seed.align
	rm Silva.seed_v123.tgz silva.seed_v123.*

$(REFS)/silva.v4.align : $(REFS)/silva.seed.align
	mothur "#pcr.seqs(fasta=$(REFS)/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
	mv $(REFS)/silva.seed.pcr.align $(REFS)/silva.v4.align

# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

$(REFS)/trainset14_032015.% :
	wget -N http://mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
	tar xvzf Trainset14_032015.pds.tgz trainset14_032015.pds/trainset14_032015.pds.*
	mv trainset14_032015.pds/* $(REFS)/
	rmdir trainset14_032015.pds
	rm Trainset14_032015.pds.tgz

################################################################################
#
# Part 2: Run data through mothur
#
#	Process fastq data through the generation of files that will be used in the
# overall analysis.
#
################################################################################

# Change gf_cdiff to the * part of your *.files file that lives in data/raw/
STUB = $(foreach S, $(STUDIES), $(PROC)/$(S)/$(S))


ALPHA = $(addsuffix .groups.ave-std.summary,$(STUB))
BETA = $(addsuffix .braycurtis.0.03.lt.ave.dist,$(STUB))
SHARED = $(addsuffix .shared,$(STUB))
SUBSHARED = $(addsuffix .0.03.subsample.shared,$(STUB))
FASTA = $(addsuffix .rep.seqs,$(STUB))
TAXONOMY = $(addsuffix .taxonomy,$(STUB))
METADATA = $(addsuffix .metadata,$(STUB))


.SECONDEXPANSION:
data/process/%.groups.ave-std.summary\
	data/process/%.braycurtis.0.03.lt.ave.dist\
	data/process/%.shared\
	data/process/%.0.03.subsample.shared\
	data/process/%.rep.seqs\
	data/process/%.taxonomy\
	data/process/%.metadata : code/$$(notdir $$*).batch code/$$(notdir $$*).R\
			$(REFS)/silva.seed.align $(REFS)/silva.v4.align\
			$(REFS)/trainset14_032015.pds.fasta\
			$(REFS)/trainset14_032015.pds.tax
	bash $<




################################################################################
#
# Part 3: Metadata Processing and Analysis
#
#	Run scripts that analyze the generated data
#
################################################################################

# Set up output files from power transformation
COMMON_PATH=$(addprefix $(TABLES)/,$(STUDIES))
ALPHA_STOOL_RAW=$(addsuffix _stool_alpha_raw_values.csv,$(COMMON_PATH))
ALPHA_STOOL_SUMMARY=$(addsuffix _stool_summary_stats_alpha_raw_values.csv,$(COMMON_PATH))
ALPHA_STOOL_TRANS=$(addsuffix _stool_transformed_data_alpha_raw_values.csv,$(COMMON_PATH))
ALPHA_TISS_RAW=$(addsuffix _tissue_alpha_raw_values.csv,$(COMMON_PATH))
ALPHA_TISS_SUMMARY=$(addsuffix _tissue_summary_stats_alpha_raw_values.csv,$(COMMON_PATH))
ALPHA_TISS_TRANS=$(addsuffix _tissue_transformed_data_alpha_raw_values.csv,$(COMMON_PATH))

# Code to power transform the Alpha diversity values of each study
$(ALPHA_STOOL_RAW)\
$(ALPHA_STOOL_SUMMARY)\
$(ALPHA_STOOL_TRANS)\
$(ALPHA_TISS_RAW)\
$(ALPHA_TISS_SUMMARY)\
$(ALPHA_TISS_TRANS)\
$(TABLES)/stool_power_transformation_summary.csv : $(METADATA) $(ALPHA)\
code/get_power_transform.R
	R -e "source('code/get_power_transform.R')"


# Code to run the stool Alpha analysis
$(TABLES)/alpha_tukey_results.csv\
$(TABLES)/alpha_ttest_results.csv\
$(TABLES)/alpha_mixeffect_results.csv\
$(TABLES)/stool_normalized_alpha_all.csv : $(ALPHA_STOOL_TRANS)\
code/get_stool_combined_alpha.R
	R -e "source('code/get_stool_combined_alpha.R')"


#Generate matched and unmatched tissue samples
$(TABLES)/alpha_tissue_matched_data.csv\
$(TABLES)/alpha_tissue_unmatched_data.csv : $(METADATA) $(ALPHA_TISS_TRANS)\
$(ALPHA) $(PROC)/geng/GengData.txt $(PROC)/dejea/SraRunTable.txt\
$(PROC)/kostic/SraRunTable.txt $(PROC)/burns/burnsMetadata.csv\
$(PROC)/lu/luData.csv code/get_matched_tissue.R
	R -e "source('code/get_matched_tissue.R')"

# Code to run the tissue Alpha analysis
$(TABLES)/alpha_unmatched_ttest_tissue.csv\
$(TABLES)/alpha_matched_ttest_tissue.csv\
$(TABLES)/alpha_unmatched_tukey_results.csv\
$(TABLES)/alpha_unmatched_mixeffect_results.csv\
$(TABLES)/alpha_matched_mixeffect_results.csv\
$(TABLES)/matched_tissue_normalized_alpha_all_data.csv\
$(TABLES)/unmatched_tissue_normalized_alpha_all_data.csv : $(TABLES)/alpha_tissue_matched_data.csv\
$(TABLES)/alpha_tissue_unmatched_data.csv code/get_tissue_combined_alpha.R
	R -e "source('code/get_tissue_combined_alpha.R')"


# Code to run stool Beta analysis - adneoma and carcinoma
$(TABLES)/beta_perm_stool_summary.csv\
$(TABLES)/beta_perm_adn_stool_summary.csv : $(BETA) $(ALPHA_STOOL_TRANS)\
code/get_beta_stool.R code/get_adenoma_beta_stool.R
	R -e "source('code/get_beta_stool.R')"
	R -e "source('code/get_adenoma_beta_stool.R')"

# Code to run the tissue Beta analysis - adenoma and carcinoma
$(TABLES)/beta_perm_unmatched_tissue_summary.csv\
$(TABLES)/beta_perm_matched_tissue_summary.csv\
$(TABLES)/bray_matched_v_cont_tissue_summary.csv\
$(TABLES)/bray_matched_v_cases_tissue_summary.csv\
$(TABLES)/beta_perm_adn_unmatched_tissue_summary.csv\
$(TABLES)/beta_perm_adn_matched_tissue_summary.csv\
$(TABLES)/bray_adn_matched_v_cont_tissue_summary.csv\
$(TABLES)/bray_adn_matched_v_cases_tissue_summary.csv : $(TABLES)/alpha_tissue_matched_data.csv\
$(TABLES)/alpha_tissue_unmatched_data.csv $(BETA)\
code/get_beta_tissue.R code/get_adenoma_beta_tissue.R
	R -e "source('code/get_beta_tissue.R')"
	R -e "source('code/get_adenoma_beta_tissue.R')"


# Generate Relative Risk for stool alpha diversity
$(TABLES)/alpha_group_counts_summary.csv\
$(TABLES)/alpha_OR_ind_results.csv\
$(TABLES)/alpha_OR_composite.csv\
$(TABLES)/alpha_adn_group_counts_summary.csv\
$(TABLES)/alpha_adn_OR_ind_results.csv\
$(TABLES)/alpha_adn_OR_composite.csv : $(ALPHA_STOOL_RAW)\
code/get_stool_RRisk.R code/get_adenoma_stool_RRisk.R
	R -e "source('code/get_stool_OR.R')"
	R -e "source('code/get_adenoma_stool_OR.R')"


# Generate Relative Risk for tissue alpha diversity
$(TABLES)/alpha_group_counts_unmatched_tissue_summary.csv\
$(TABLES)/alpha_OR_ind_unmatched_tissue_results.csv\
$(TABLES)/alpha_OR_unmatched_tissue_composite.csv\
$(TABLES)/alpha_group_counts_tissue_summary.csv\
$(TABLES)/alpha_OR_ind_tissue_results.csv\
$(TABLES)/alpha_OR_tissue_composite.csv\
$(TABLES)/alpha_group_counts_matched_tissue_summary.csv\
$(TABLES)/alpha_OR_ind_matched_tissue_results.csv\
$(TABLES)/alpha_OR_matched_tissue_composite.csv\
$(TABLES)/alpha_adn_group_counts_tissue_summary.csv\
$(TABLES)/alpha_adn_OR_ind_tissue_results.csv\
$(TABLES)/alpha_adn_OR_tissue_composite.csv : $(TABLES)/alpha_tissue_matched_data.csv\
$(TABLES)/alpha_tissue_unmatched_data.csv code/get_tissue_RRisk.R\
code/get_adenoma_tissue_RRisk.R
	R -e "source('code/get_tissue_OR.R')"
	R -e "source('code/get_adenoma_tissue_OR.R')"


# Set up the genera share file variable
GENERA_FILE=$(addsuffix _genera_shared.csv,$(STUB))
SUB_GENERA_FILE=$(addsuffix _subsample_genera.csv,$(STUB))

# Create Genus level shared file
$(GENERA_FILE) : $(ALPHA) $(TAXONOMY) $(SHARED)\
code/get_taxonomy.R
	R -e "source('code/get_taxonomy.R')"

# Create a subsampled Genus level shared file
$(SUB_GENERA_FILE) : $(GENERA_FILE) code/get_genera_randomization.R
	R -e "source('code/get_genera_randomization.R')"


# Analyze Stool Genera Relative Risk
$(TABLES)/select_genus_group_counts_summary.csv\
$(TABLES)/select_genus_OR_stool_ind_results.csv\
$(TABLES)/select_genus_OR_stool_composite.csv\
$(TABLES)/adn_select_genus_group_counts_summary.csv\
$(TABLES)/adn_select_genus_OR_stool_ind_results.csv\
$(TABLES)/adn_select_genus_OR_stool_composite.csv : $(GENERA_FILE) $(SUB_GENERA_FILE)\
$(METADATA) code/get_specifc_genus_analysis_stool.R\
code/get_select_inc_genera_positivity_OR_stool.R\
code/get_adenoma_specifc_genus_analysis_stool.R\
code/get_adneoma_select_inc_genera_positivity_OR_stool.R
	R -e "source('code/get_specifc_genus_analysis_stool.R')"
	R -e "source('code/get_adenoma_specifc_genus_analysis_stool.R')"



# Analyze Tissue Genera Relative Risk
$(TABLES)/select_genus_matched_tissue_group_counts_summary.csv\
$(TABLES)/select_genus_OR_matched_tissue_ind_results.csv\
$(TABLES)/select_genus_OR_matched_tissue_composite.csv\
$(TABLES)/select_genus_unmatched_tissue_group_counts_summary.csv\
$(TABLES)/select_genus_OR_unmatched_tissue_ind_results.csv\
$(TABLES)/select_genus_OR_unmatched_tissue_composite.csv\
$(TABLES)/select_genus_tissue_group_counts_summary.csv\
$(TABLES)/select_genus_OR_tissue_ind_results.csv\
$(TABLES)/select_genus_OR_tissue_composite.csv\
$(TABLES)/adn_select_genus_tissue_group_counts_summary.csv\
$(TABLES)/adn_select_genus_OR_tissue_ind_results.csv\
$(TABLES)/adn_select_genus_OR_tissue_composite : $(GENERA_FILE) $(SUB_GENERA_FILE)\
$(TABLES)/alpha_tissue_matched_data.csv $(TABLES)/alpha_tissue_unmatched_data.csv\
code/get_specifc_genus_analysis_tissue.R\
code/get_select_inc_genera_positivity_OR_tissue.R\
code/get_adenoma_specific_genus_analysis_tissue.R\
code/get_adneoma_select_inc_genera_positivity_OR_tissue.R
	R -e "source('code/get_specifc_genus_analysis_tissue.R')"
	R -e "source('code/get_adenoma_specific_genus_analysis_tissue.R')"


# Create Variables for the stool dependencies
CRC_STOOL_STUDY = wang weir ahn zeller baxter hale flemer
ADN_STOOL_STUDY = brim zeller baxter hale

# Set up taxonomy call for important variable counting
CRC_STOOL_TAX = $(foreach S, $(CRC_STOOL_STUDY), $(PROC)/$(S).taxonomy)
ADN_STOOL_TAX = $(foreach S, $(ADN_STOOL_STUDY), $(PROC)/$(S).taxonomy)

# Set up directory for stool genera RF
G_CRC_STOOL_FULL = $(foreach S, $(CRC_STOOL_STUDY), $(TABLES)/genus_stool_RF_full_$(S))
G_CRC_STOOL_SELECT = $(foreach S, $(CRC_STOOL_STUDY), $(TABLES)/genus_stool_RF_select_$(S))
G_ADN_STOOL_FULL = $(foreach S, $(ADN_STOOL_STUDY), $(TABLES)/adn_genus_stool_RF_full_$(S))
G_ADN_STOOL_SELECT = $(foreach S, $(ADN_STOOL_STUDY), $(TABLES)/adn_genus_stool_RF_select_$(S))

# Set up files to be created for stool genera RF
G_CRC_FULL_STOOL_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_CRC_STOOL_FULL))
G_CRC_FULL_STOOL_ROC=$(addsuffix _raw_roc_data.csv,$(G_CRC_STOOL_FULL))
G_CRC_FULL_STOOL_IMP=$(addsuffix _imp_vars.csv,$(G_CRC_STOOL_FULL))
G_CRC_SELECT_STOOL_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_CRC_STOOL_SELECT))
G_CRC_SELECT_STOOL_ROC=$(addsuffix _raw_roc_data.csv,$(G_CRC_STOOL_SELECT))

G_ADN_FULL_STOOL_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_ADN_STOOL_FULL))
G_ADN_FULL_STOOL_ROC=$(addsuffix _raw_roc_data.csv,$(G_ADN_STOOL_FULL))
G_ADN_FULL_STOOL_IMP=$(addsuffix _imp_vars.csv,$(G_ADN_STOOL_FULL))
G_ADN_SELECT_STOOL_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_ADN_STOOL_SELECT))
G_ADN_SELECT_STOOL_ROC=$(addsuffix _raw_roc_data.csv,$(G_ADN_STOOL_SELECT))

# Set up files to be created for OTU RF
O_CRC_STOOL_FULL = $(foreach S, $(CRC_STOOL_STUDY), $(TABLES)/$(S))
O_ADN_STOOL_FULL = $(foreach S, $(ADN_STOOL_STUDY), $(TABLES)/adn_$(S))

O_CRC_FULL_STOOL_IMP=$(addsuffix _imp_otu_table.csv,$(O_CRC_STOOL_FULL))
O_ADN_FULL_STOOL_IMP=$(addsuffix _imp_otu_table.csv,$(O_ADN_STOOL_FULL))

# Run the Stool Genera Random Forest Models
$(G_CRC_FULL_STOOL_PVALUE) $(G_CRC_FULL_STOOL_ROC)\
$(G_CRC_SELECT_STOOL_PVALUE) $(G_CRC_SELECT_STOOL_ROC)\
$(G_CRC_FULL_STOOL_IMP) $(G_ADN_FULL_STOOL_IMP)\
$(G_ADN_FULL_STOOL_PVALUE) $(G_ADN_FULL_STOOL_ROC)\
$(G_ADN_SELECT_STOOL_PVALUE) $(G_ADN_SELECT_STOOL_ROC)\
$(TABLES)/genus_stool_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/adn_genus_stool_RF_fullvsselect_pvalue_summary.csv : $(GENERA_FILE)\
$(SUB_GENERA_FILE) $(METADATA)\
code/run_random_forest_genus_stool.R\
code/run_adn_random_forest_genus_stool.R
	R -e "source('code/run_non_common_RF_genera_stool.R')"
	R -e "source('code/run_adn_non_common_RF_genera_stool.R')"


# Create Variables for the tissue dependencies
CRC_MATCH_TISSUE_STUDY = burns dejea geng 
CRC_UNMATCH_TISSUE_STUDY = burns chen flemer 
ADN_TISSUE_STUDY = lu flemer sana

# Set up taxonomy call for important variable counting
CRC_MATCH_T_TAX = $(foreach S, $(CRC_MATCH_TISSUE_STUDY), $(PROC)/$(S).taxonomy)
CRC_UNMATCH_T_TAX = $(foreach S, $(CRC_UNMATCH_TISSUE_STUDY), $(PROC)/$(S).taxonomy)
ADN_T_TAX = $(foreach S, $(ADN_TISSUE_STUDY), $(PROC)/$(S).taxonomy)

# Set up directory for genera tissue RF
G_CRC_MATCH_T_FULL = $(foreach S, $(CRC_MATCH_TISSUE_STUDY), $(TABLES)/genus_matched_tissue_RF_$(S))
G_CRC_MATCH_T_SELECT = $(foreach S, $(CRC_MATCH_TISSUE_STUDY), $(TABLES)/genus_matched_tissue_RF_select_$(S))

G_CRC_UNMATCH_T_FULL = $(foreach S, $(CRC_UNMATCH_TISSUE_STUDY), $(TABLES)/genus_unmatched_tissue_RF_$(S))
G_CRC_UNMATCH_T_SELECT = $(foreach S, $(CRC_UNMATCH_TISSUE_STUDY), $(TABLES)/genus_unmatched_tissue_RF_select_$(S))

G_ADN_TISSUE_FULL = $(foreach S, $(ADN_TISSUE_STUDY), $(TABLES)/adn_genus_unmatched_tissue_RF_full_$(S))
G_ADN_TISSUE_SELECT = $(foreach S, $(ADN_TISSUE_STUDY), $(TABLES)/adn_genus_unmatched_tissue_RF_select_$(S))

# Set up files to be created for genera tissue RF
G_CRC_FULL_MATCH_T_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_CRC_MATCH_T_FULL))
G_CRC_FULL_MATCH_T_ROC=$(addsuffix _raw_roc_data.csv,$(G_CRC_MATCH_T_FULL))
G_CRC_FULL_MATCH_T_IMP=$(addsuffix _imp_vars.csv,$(G_CRC_MATCH_T_FULL))
G_CRC_FULL_UNMATCH_T_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_CRC_UNMATCH_T_FULL))
G_CRC_FULL_UNMATCH_T_ROC=$(addsuffix _raw_roc_data.csv,$(G_CRC_UNMATCH_T_FULL))
G_CRC_FULL_UNMATCH_T_IMP=$(addsuffix _imp_vars.csv,$(G_CRC_UNMATCH_T_FULL))

G_CRC_SELECT_MATCH_T_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_CRC_MATCH_T_SELECT))
G_CRC_SELECT_MATCH_T_ROC=$(addsuffix _raw_roc_data.csv,$(G_CRC_MATCH_T_SELECT))
G_CRC_SELECT_UNMATCH_T_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_CRC_UNMATCH_T_SELECT))
G_CRC_SELECT_UNMATCH_T_ROC=$(addsuffix _raw_roc_data.csv,$(G_CRC_UNMATCH_T_SELECT))


G_ADN_FULL_TISSUE_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_ADN_TISSUE_FULL))
G_ADN_FULL_TISSUE_ROC=$(addsuffix _raw_roc_data.csv,$(G_ADN_TISSUE_FULL))
G_ADN_FULL_TISSUE_IMP=$(addsuffix _imp_vars.csv,$(G_ADN_TISSUE_FULL))
G_ADN_SELECT_TISSUE_PVALUE=$(addsuffix _pvalue_summary.csv,$(G_ADN_TISSUE_SELECT))
G_ADN_SELECT_TISSUE_ROC=$(addsuffix _raw_roc_data.csv,$(G_ADN_TISSUE_SELECT))

# Set up files to be created for OTU RF
O_CRC_MATCH_T_FULL = $(foreach S, $(CRC_MATCH_TISSUE_STUDY), $(TABLES)/$(S))
O_CRC_UNMATCH_T_FULL = $(foreach S, $(CRC_UNMATCH_TISSUE_STUDY), $(TABLES)/$(S))
O_ADN_TISSUE_FULL = $(foreach S, $(ADN_TISSUE_STUDY), $(TABLES)/adn_$(S))

O_CRC_FULL_MATCH_T_IMP=$(addsuffix _matched_tissue_imp_otu_table.csv,$(O_CRC_MATCH_T_FULL))
O_CRC_FULL_UNMATCH_T_IMP=$(addsuffix _unmatched_tissue_imp_otu_table.csv,$(O_CRC_UNMATCH_T_FULL))
O_ADN_TISSUE_FULL_T_IMP=$(addsuffix _tissue_imp_otu_table.csv,$(O_ADN_TISSUE_FULL))

# Run the Tissue Genera Random Forest Models
$(G_CRC_FULL_MATCH_T_PVALUE) $(G_CRC_FULL_MATCH_T_ROC)\
$(G_CRC_FULL_UNMATCH_T_PVALUE) $(G_CRC_FULL_UNMATCH_T_ROC)\
$(G_CRC_SELECT_MATCH_T_PVALUE) $(G_CRC_SELECT_MATCH_T_ROC)\
$(G_CRC_SELECT_UNMATCH_T_PVALUE) $(G_CRC_SELECT_UNMATCH_T_ROC)\
$(G_ADN_FULL_TISSUE_PVALUE) $(G_ADN_FULL_TISSUE_ROC)\
$(G_ADN_SELECT_TISSUE_PVALUE) $(G_ADN_SELECT_TISSUE_ROC)\
$(G_CRC_FULL_MATCH_T_IMP) $(G_CRC_FULL_UNMATCH_T_IMP)\
$(G_ADN_FULL_TISSUE_IMP)\
$(TABLES)/genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/genus_matched_tissue_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/adn_genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)adn_genus_matched_tissue_RF_fullvsselect_pvalue_summary.csv : $(GENERA_FILE)\
$(SUB_GENERA_FILE) $(TABLES)/alpha_tissue_matched_data.csv\
$(TABLES)/alpha_tissue_unmatched_data.csv\
code/run_random_forest_genus_tissue.R code/run_adn_random_forest_genus_tissue.R
	R -e "source('code/run_non_common_RF_genera_tissue.R')"
	R -e "source('code/run_adn_non_common_RF_genera_tissue.R')"



#$(O_ADN_FULL_STOOL_IMP)\

# Run the Stool OTU Random Forest Models
$(TABLES)/stool_rf_otu_roc.csv\
$(TABLES)/stool_rf_otu_random_comparison_summary.csv\
$(TABLES)/stool_rf_select_otu_roc.csv\
$(TABLES)/stool_rf_select_otu_random_comparison_summary.csv\
$(TABLES)/adn_stool_rf_otu_roc.csv\
$(O_CRC_FULL_STOOL_IMP) $(O_ADN_FULL_STOOL_IMP)\
$(TABLES)/adn_stool_rf_otu_random_comparison_summary.csv : $(SUBSHARED) $(METADATA)\
code/run_random_forest_otu_stool.R\
code/run_adn_random_forest_otu_stool.R\
code/run_random_forest_select_otu_stool.R\
code/run_adn_random_forest_select_otu_stool.R
	R -e "source('code/run_random_forest_otu_stool.R')"
	R -e "source('code/run_random_forest_select_otu_stool.R')"
	R -e "source('code/run_adn_random_forest_otu_stool.R')"
	R -e "source('code/run_adn_random_forest_select_otu_stool.R')"



# Run the Tissue OTU Random Forest Models
$(TABLES)/unmatched_tissue_rf_otu_roc.csv\
$(TABLES)/unmatched_tissue_rf_otu_random_comparison_summary.csv\
$(TABLES)/matched_tissue_rf_otu_roc.csv\
$(TABLES)/matched_tissue_rf_otu_random_comparison_summary.csv\
$(TABLES)/adn_tissue_rf_otu_roc.csv\
$(O_CRC_FULL_MATCH_T_IMP) $(O_CRC_FULL_UNMATCH_T_IMP)\
$(O_ADN_TISSUE_FULL_T_IMP)\
$(TABLES)/adn_tissue_rf_otu_random_comparison_summary.csv : $(SUBSHARED) S(METADATA)\
$(TABLES)/alpha_tissue_matched_data.csv $(TABLES)/alpha_tissue_unmatched_data.csv\
code/run_random_forest_otu_tissue.R code/run_adn_random_forest_otu_tissue.R\
code/run_random_forest_select_otus_tissue.R code/run_adn_random_forest_select_otu_tissue
	R -e "source('code/run_random_forest_otu_tissue.R')"
	R -e "source('code/run_random_forest_select_otus_tissue.R')"
	R -e "source('code/run_adn_random_forest_otu_tissue.R')"
	R -e "source('code/run_adn_random_forest_select_otu_tissue.R')"



# Generate the most important variables with the RF across study for stool
$(TABLES)/crc_RF_genera_stool_top10.csv\
$(TABLES)/adn_RF_genera_stool_top10.csv\
$(TABLES)/crc_RF_otu_stool_top10.csv\
$(TABLES)/adn_RF_otu_stool_top10.csv : $(CRC_STOOL_TAX) $(ADN_STOOL_TAX)\
$(G_CRC_FULL_STOOL_IMP) $(G_ADN_FULL_STOOL_IMP) $(O_CRC_FULL_STOOL_IMP)\
$(O_ADN_FULL_STOOL_IMP) code/run_comparison_RF_imp_otus_stool.R
	R -e "source('code/run_comparison_RF_imp_otus_stool.R')"


# Generate the most important variables with the RF across study for tissue
$(TABLES)/crc_RF_genera_unmatched_tissue_top10.csv\
$(TABLES)/adn_RF_genera_matched_tissue_top10.csv\
$(TABLES)/adn_RF_genera_tissue_top10.csv\
$(TABLES)/crc_RF_otu_unmatched_tissue_top10.csv\
$(TABLES)/adn_RF_otu_matched_tissue_top10.csv\
$(TABLES)/adn_RF_otu_tissue_top10.csv : $(TABLES)/alpha_tissue_matched_data.csv\
$(TABLES)/alpha_tissue_unmatched_data.csv $(CRC_MATCH_T_TAX) $(CRC_UNMATCH_T_TAX)\
$(ADN_T_TAX) $(G_CRC_FULL_MATCH_T_IMP) $(G_CRC_FULL_UNMATCH_T_IMP) $(G_ADN_FULL_TISSUE_IMP)\
$(O_CRC_FULL_MATCH_T_IMP) $(O_CRC_FULL_UNMATCH_T_IMP) $(O_ADN_FULL_STOOL_IMP)\
code/run_comparison_RF_imp_otus_tissue.R
	R -e "source('code/run_comparison_RF_imp_otus_tissue.R')"

# Generate individual taxa AUCs
$(TABLES)/ind_genera_auc_stool.csv\
$(TABLES)/ind_genera_auc_unmatched_tissue.csv : $(GENERA_FILE) $(SUB_GENERA_FILE)\
$(METADATA) $(TABLES)/select_genus_OR_stool_composite.csv $(TABLES)/alpha_tissue_unmatched_data.csv\
$(TABLES)/select_genus_OR_unmatched_tissue_composite.csv code/Run_sig_taxa_AUC_generation.R
	R -e "source('code/Run_sig_taxa_AUC_generation.R')"


# Run the pvalue AUC analysis between ind taxa and full taxa RF models
$(TABLES)/stool_ind_vs_full_taxa_pvalue_summary.csv\
$(TABLES)/unmatched_tissue_ind_vs_full_taxa_pvalue_summary.csv : $(TABLES)/ind_genera_auc_stool.csv\
$(TABLES)/ind_genera_auc_unmatched_tissue.csv\
$(TABLES)/ALL_genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/ALL_genus_matched_tissue_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/matched_tissue_rf_otu_random_comparison_summary.csv\
$(TABLES)/unmatched_tissue_rf_otu_random_comparison_summary.csv\
$(TABLES)/ALL_genus_stool_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/stool_rf_otu_random_comparison_summary.csv\
code/run_pvalue_auc_ind_vs_full_taxa.R
	R -e "source('code/run_pvalue_auc_ind_vs_full_taxa.R')"


# Generate total n used in study for stool and tissue
$(TABLES)/stool_study_n_analyzed.csv\
$(TABLES)/tissue_study_n_analyzed.csv : $(TABLES)/alpha_adn_group_counts_summary.csv\
$(TABLES)/alpha_group_counts_summary.csv $(TABLES)/alpha_adn_group_counts_tissue_summary.csv\
$(TABLES)/alpha_group_counts_matched_tissue_summary.csv $(TABLES)/alpha_group_counts_unmatched_tissue_summary.csv\
code/make_study_count_table.R
	R -e "source('code/make_study_count_table.R')"



################################################################################
#
# Part 3: Figure and table generation
#
#	Run scripts to generate figures and tables
#
################################################################################


# Run Alpha raw graph Figure 1
$(FIGS)/Figure1.pdf : $(ALPHA_STOOL_RAW) $(ALPHA_TISS_RAW)\
code/make_stool_alpha_raw_graph.R
	R -e "source('code/make_stool_alpha_raw_graph.R')"


# Run code to create Figure 2
$(FIGS)/Figure2.pdf : $(TABLES)/alpha_adn_RR_composite.csv\
$(TABLES)/alpha_adn_RR_ind_results.csv $(TABLES)/alpha_RR_composite.csv\
$(TABLES)/alpha_RR_ind_results.csv code/make_stool_alpha_RR_graph.R
	R -e "source('code/make_stool_alpha_OR_graph.R')"


# Run code to create Figure 3
$(FIGS)/Figure3.pdf : $(TABLES)/ind_genera_auc_stool.csv\
$(TABLES)/ind_genera_auc_unmatched_tissue.csv\
code/make_ind_genera_auc_graph.R
	R -e "source('code/make_ind_genera_auc_graph.R')"


# Run code to create Figure 4
$(FIGS)/Figure4.pdf : $(TABLES)/ALL_genus_stool_RF_select_imp_vars.csv\
$(TABLES)/ALL_genus_unmatched_tissue_RF_select_imp_vars.csv code/make_select_imp_otu_graph.R
	R -e "source('code/make_select_imp_otu_graph.R')"


# Run code to create Figure 5 and S2
$(FIGS)/FigureS2.pdf\
$(FIGS)/Figure5.pdf : $(TABLES)/adn_ALL_genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/adn_tissue_rf_otu_random_comparison_summary.csv\
$(TABLES)/adn_ALL_genus_stool_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/adn_stool_rf_otu_random_comparison_summary.csv\
$(TABLES)/ALL_genus_matched_tissue_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/ALL_genus_unmatched_tissue_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/matched_tissue_rf_otu_random_comparison_summary.csv\
$(TABLES)/unmatched_tissue_rf_otu_random_comparison_summary.csv\
$(TABLES)/ALL_genus_stool_RF_fullvsselect_pvalue_summary.csv\
$(TABLES)/stool_rf_otu_random_comparison_summary.csv\
code/make_rf_auc_full_versus_specific_graph.R
	R -e "source('code/make_rf_auc_full_versus_specific_graph.R')"



# Run code to make supplemental Figure S1
$(FIGS)/FigureS1.pdf : $(TABLES)/alpha_adn_OR_tissue_composite.csv\
$(TABLES)/alpha_adn_OR_ind_tissue_results.csv $(TABLES)/alpha_OR_tissue_composite.csv\
$(TABLES)/alpha_OR_ind_tissue_results.csv code/make_tissue_alpha_OR_graph.R
	R -e "source('code/make_tissue_alpha_OR_graph.R')"



# Run code to create Figure S3
$(FIGS)/FigureS3.pdf : $(TABLES)/crc_RF_genera_stool_top10.csv\
$(TABLES)/adn_RF_genera_stool_top10.csv\
$(TABLES)/crc_RF_otu_stool_top10.csv\
$(TABLES)/adn_RF_otu_stool_top10.csv code/make_stool_imp_otu_graph.R
	R -e "source('code/make_stool_imp_otu_graph.R')"

# Run code to make supplemental Figure S4
$(FIGS)/FigureS4.pdf : $(TABLES)/crc_RF_genera_unmatched_tissue_top10.csv\
$(TABLES)/adn_RF_genera_matched_tissue_top10.csv\
$(TABLES)/adn_RF_genera_tissue_top10.csv\
$(TABLES)/crc_RF_otu_unmatched_tissue_top10.csv\
$(TABLES)/adn_RF_otu_matched_tissue_top10.csv\
$(TABLES)/adn_RF_otu_tissue_top10.csv code/make_tissue_imp_otu_graph.R
	R -e "source('code/make_tissue_imp_otu_graph.R')"


# Run code to make Figure 6, S5
$(FIGS)/Figure6.pdf\
$(FIGS)/FigureS5.pdf : $(G_ADN_FULL_STOOL_PVALUE)\
$(G_ADN_FULL_TISSUE_PVALUE) $(G_CRC_FULL_STOOL_PVALUE)\
$(G_CRC_SELECT_STOOL_PVALUE) $(G_CRC_FULL_MATCH_T_PVALUE)\
$(G_CRC_FULL_UNMATCH_T_PVALUE) $(G_CRC_SELECT_UNMATCH_T_PVALUE)\
code/make_genus_rf_auc_against_study_graph.R
	R -e "source('code/make_genus_rf_auc_against_study_graph.R')"



################################################################################
#
# Part 4: Pull it all together
#
# Render the manuscript
#
################################################################################

write.paper : $(FINAL)/manuscript.Rmd $(FINAL)/supplement.Rmd\
$(FIGS)/Figure1.pdf $(FIGS)/Figure2.pdf\
$(FIGS)/Figure3.pdf $(FIGS)/Figure4.pdf\
$(FIGS)/Figure5.pdf $(FIGS)/FigureS1.pdf\
$(FIGS)/FigureS2.pdf $(FIGS)/FigureS3.pdf\
$(FIGS)/FigureS4.pdf code/Run_render_paper.R
	R -e "source('code/Run_render_paper.R')"


submission/response_to_reviewers.pdf : submission/response_to_reviewers.md
	pandoc -s --include-in-header=submission/header.tex -V geometry:margin=1in -o $@ $<
	
	
write.revision1.paper : $(FINAL)/manuscript_R1.Rmd\
$(FIGS)/Figure1.pdf $(FIGS)/Figure2.pdf\
$(FIGS)/Figure3.pdf $(FIGS)/Figure4.pdf\
$(FIGS)/Figure5.pdf $(FIGS)/FigureS1.pdf\
$(FIGS)/FigureS2.pdf $(FIGS)/FigureS3.pdf\
$(FIGS)/FigureS4.pdf code/Run_render_paper.R
	R -e "source('code/Run_render_revision1_paper.R')"

write.r1.marked.up : $(FINAL)/manuscript.tex\
$(FINAL)/manuscript_R1.tex
	latexdiff $(FINAL)/manuscript.tex $(FINAL)/manuscript_R1.tex > $(FINAL)/manuscript_R1_markedup.tex
	pdflatex -output-directory=$(FINAL) $(FINAL)/manuscript_R1_markedup.tex
	rm $(FINAL)/manuscript_R1_markedup.aux 
	rm $(FINAL)/manuscript_R1_markedup.log
	rm $(FINAL)/manuscript_R1_markedup.out	

