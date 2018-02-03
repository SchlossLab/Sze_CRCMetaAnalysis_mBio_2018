## Making Sense of the Noise: Leveraging Existing 16S rRNA Gene Surveys to Identify Key Community Members in Colorectal Tumors


**Background.** An increasing body of literature suggests that there is a crucial role for the microbiota in colorectal cancer (CRC) pathogenesis. Important drivers within this context have ranged from individual microbes to the whole community. Our study expands on a recent meta-analysis investigating microbial biomarkers for CRC by testing the hypothesis that the bacterial community has important associations to both early (adenoma) and late (carcinoma) stage disease. To test this hypothesis we examined both feces (n = 1737) and colon tissue (492 total samples from 350 individuals) across 14 previously published 16S rRNA gene sequencing studies on CRC and the microbiota.


**Results.** Fecal samples had a significant decrease for both Shannon diversity and
evenness, after correcting for study effect and variable region sequenced, with more
severe disease (P-value < 0.05). This reduction in evenness translated into small increases in relative risk for adenoma (P-value = 0.032) and carcinoma stages of CRC (P-value = 0.00034) while the reduction in Shannon diversity only translated into an increased relative risk for developing carcinomas (P-value = 0.0047). Increases in mouth-associated microbes were commonly in the top 5 most significantly increased relative risk of adenoma and carcinoma for both stool and tissue samples. A prediction model for adenoma and carcinoma was built using either the whole community or selected genera with highest and lowest relative risk from fecal and tissue samples. Both approaches resulted in similar classification success according to Area Under the Curve (AUC) regardless of whether genera or OTUs were used to build the model. The most important groups within the full community models consistently belonged to genera such as Ruminococcus, Bacteroides, and Roseburia across studies. Although a number of associations between the microbiota and CRC were identified, the majority of studies that we used in this meta-analysis were
only individually adequately powered for large effect sizes.


**Conclusions.** These data provide support for the importance of the bacterial community to both adenoma and carcinoma genesis. The evidence collected within this study on the role of the microbiota in CRC identifies a number of correlations that may not have been detected because of the low power associated with the majority of studies that have been performed to date.




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
