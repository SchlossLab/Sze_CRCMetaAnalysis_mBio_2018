## Investigating the Microbiota and Colorectal Cancer: The Importance of Community

**Background.** An increasing body of literature suggests that there is a crucial role for the
microbiota in colorectal cancer (CRC) pathogenesis. Important drivers within this context
have ranged from individual microbes to the whole community. Our study expands on a
recent meta-analysis investigating microbial biomarkers for CRC by testing the hypothesis
that the bacterial community is an important driver of both early (adenoma) and late
(carcinoma) stage of disease. To test this hypothesis we examined both feces (n = 1737)
and tissue (492 total samples from 350 individuals) across 14 different studies.

**Results.** Fecal samples had a significant decrease from control to adenoma to carcinoma
for both Shannon diversity and evenness (P-value < 0.05) after correcting for study effect
and variable region sequenced. Only evenness for adenoma (P-value < 0.05) resulted in a
slightly increased relative risk while lower Shannon diversity and evenness in fecal samples
resulted in a significant increase in relative risk for carcinoma (P-value < 0.05). Previously
associated colorectal cancer genera (Fusobacterium, Parvimonas, Peptostreptococcus,
or Porphyromonas) followed a similar pattern with a significantly increased relative risk
by their presence for carcinoma (P-value < 0.05) but not adenoma (P-value > 0.05) with
the exception of Porphyromonas (P-value < 0.05). Using the whole community versus
only CRC associated genera to build a prediction model resulted in higher classification
success based on Area Under the Curve (AUC) for both adenoma and carcinoma using
fecal and tissue samples. The most important OTUs for these models consistently belonged to genera such as *Ruminococcus*, *Bacteroides*, and *Roseburia* across studies. For the included studies, most were adequately powered for large effect size differences. This may be more amenable for carcinoma than for adenoma microbiota research due to the smaller community level changes observed in our results.

**Conclusions.** This data provides support for the importance of the bacteral community to
both adenoma and carcinoma genesis. The evidence collected within this study on the
role of the microbiota in CRC pathogenesis is much stronger for carcinoma then adenoma.
A strong reason for this may be in part due to the low power to detect more subtle changes
in the majority of studies that have been performed to date.



### Overview
	project
	|- README # the top level description of content (this doc) - 
	|CONTRIBUTING # instructions for how to contribute to your 
	|project - LICENSE # the license for this project
	|
	|- submission/
	| |- study.Rmd # executable Rmarkdown for this study, if 
	| |applicable - study.md # Markdown (GitHub) version of the 
	| |*.Rmd file - study.tex # TeX version of *.Rmd file - 
	| |study.pdf # PDF version of *.Rmd file - header.tex # LaTeX 
	| |header file to format pdf version of manuscript - 
	| |references.bib # BibTeX formatted references - XXXX.csl # csl 
	| |file to format references for journal XXX
	|
	|- data # raw and primary data, are not changed once created
	| |- references/ # reference files to be used in analysis - raw/ 
	| |# raw data, will not be altered
	| |- mothur/ # mothur processed data
	| +- process/ # cleaned data, will not be altered once created;
	|                 # will be committed to repo
	|
	|- code/ # any programmatic code
	|
	|- results # all output from workflows and analyses
	| |- tables/ # text version of tables to be rendered with kable 
	| |in R - figures/ # graphs, likely designated for manuscript 
	| |figures
	| +- pictures/ # diagrams, images, and other non-graph graphics
	|
	|- exploratory/ # exploratory data analysis for study
	| |- notebook/ # preliminary analyses
	| +- scratch/ # temporary files that can be safely deleted or 
	| lost
	|
	+- Makefile # executable Makefile for this study, if applicable

### How to regenerate this repository
#### Dependencies and locations  
* Gnu Make (v3.81) should be located in the user's PATH  
* mothur (v1.39.3) should be located in the user's PATH
	* Note v1.37.0 causes all sorts of headaches
* sratoolkit should be located in the user's PATH  	
* R (v. 3.4.2) should be located in the user's PATH  

#### Running analysis  
```git clone https://github.com/SchlossLab/Sze_CRCMetaAnalysis_XXXX_2016.git```  
```make write.paper```
