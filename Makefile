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

#$(REFS)/silva.seed.align :
#	wget -N http://mothur.org/w/images/1/15/Silva.seed_v123.tgz
#	tar xvzf Silva.seed_v123.tgz silva.seed_v123.align silva.seed_v123.tax
#	mothur "#get.lineage(fasta=silva.seed_v123.align, taxonomy=silva.seed_v123.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v123.pick.align, processors=8)"
#	mv silva.seed_v123.pick.align $(REFS)/silva.seed.align
#	rm Silva.seed_v123.tgz silva.seed_v123.*

#$(REFS)/silva.v4.align : $(REFS)/silva.seed.align
#	mothur "#pcr.seqs(fasta=$(REFS)/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
#	mv $(REFS)/silva.seed.pcr.align $(REFS)/silva.v4.align

# Next, we want the RDP reference taxonomy. The current version is v10 and we
# use a "special" pds version of the database files, which are described at
# http://blog.mothur.org/2014/10/28/RDP-v10-reference-files/

#$(REFS)/trainset14_032015.% :
#	wget -N http://mothur.org/w/images/8/88/Trainset14_032015.pds.tgz
#	tar xvzf Trainset14_032015.pds.tgz trainset14_032015.pds/trainset14_032015.pds.*
#	mv trainset14_032015.pds/* $(REFS)/
#	rmdir trainset14_032015.pds
#	rm Trainset14_032015.pds.tgz

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


#.SECONDEXPANSION:
#data/process/%.groups.ave-std.summary\
#	data/process/%.braycurtis.0.03.lt.ave.dist\
#	data/process/%.shared\
#	data/process/%.0.03.subsample.shared\
#	data/process/%.rep.seqs\
#	data/process/%.taxonomy\
#	data/process/%.metadata : code/$$(notdir $$*).batch code/$$(notdir $$*).R\
#			$(REFS)/silva.seed.align $(REFS)/silva.v4.align\
#			$(REFS)/trainset14_032015.pds.fasta\
#			$(REFS)/trainset14_032015.pds.tax
#	bash $<




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
$(TABLES)/alpha_RR_ind_results.csv\
$(TABLES)/alpha_RR_composite.csv\
$(TABLES)/alpha_adn_group_counts_summary.csv\
$(TABLES)/alpha_adn_RR_ind_results.csv\
$(TABLES)/alpha_adn_RR_composite.csv : $(ALPHA_STOOL_RAW)\
code/get_stool_RRisk.R code/get_adenoma_stool_RRisk.R
	R -e "source('code/get_stool_RRisk.R')"
	R -e "source('code/get_adenoma_stool_RRisk.R')"


# Generate Relative Risk for tissue alpha diversity
$(TABLES)/alpha_group_counts_unmatched_tissue_summary.csv\
$(TABLES)/alpha_RR_ind_unmatched_tissue_results.csv\
$(TABLES)/alpha_RR_unmatched_tissue_composite.csv\
$(TABLES)/alpha_group_counts_tissue_summary.csv\
$(TABLES)/alpha_RR_ind_tissue_results.csv\
$(TABLES)/alpha_RR_tissue_composite.csv\
$(TABLES)/alpha_group_counts_matched_tissue_summary.csv\
$(TABLES)/alpha_RR_ind_matched_tissue_results.csv\
$(TABLES)/alpha_RR_matched_tissue_composite.csv\
$(TABLES)/alpha_adn_group_counts_tissue_summary.csv\
$(TABLES)/alpha_adn_RR_ind_tissue_results.csv\
$(TABLES)/alpha_adn_RR_tissue_composite.csv : $(TABLES)/alpha_tissue_matched_data.csv\
$(TABLES)/alpha_tissue_unmatched_data.csv code/get_tissue_RRisk.R\
code/get_adenoma_tissue_RRisk.R
	R -e "source('code/get_tissue_RRisk.R')"
	R -e "source('code/get_adenoma_tissue_RRisk.R')"





################################################################################
#
# Part 3: Figure and table generation
#
#	Run scripts to generate figures and tables
#
################################################################################


################################################################################
#
# Part 4: Pull it all together
#
# Render the manuscript
#
################################################################################


#$(FINAL)/study.% : 			\ #include data files that are needed for paper
#						$(FINAL)/peerj.csl\
#						$(FINAL)/references.bib\
#						$(FINAL)/study.Rmd
#	R -e 'render("$(FINAL)/study.Rmd", clean=FALSE)'
#	mv $(FINAL)/study.knit.md $@
#	rm $(FINAL)/study.utf8.md

#write.paper : $(TABLES)/table_1.pdf $(TABLES)/table_2.pdf\ #customize to include
#				$(FIGS)/figure_1.pdf $(FIGS)/figure_2.pdf\	# appropriate tables and
#				$(FIGS)/figure_3.pdf $(FIGS)/figure_4.pdf\	# figures
#				$(FINAL)/study.Rmd $(FINAL)/study.md\
#				$(FINAL)/study.tex $(FINAL)/study.pdf
