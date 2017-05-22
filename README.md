## Colorectal Cancer Meta-Analysis of Stool and Tissue

### Main Questions to Answer

* Can we identify global discriminative OTUs in stool that classify tissue?
* Can we identify global discriminative OTUs in tissue that classify stool?
* Using all data is there a common signature of health and disease



### Active To Do list

* Need to run the Lu data set from the cluster split step onward to incorporate opticlust


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
* mothur (v1.37.4) should be located in the user's PATH
	* Note v1.37.0 causes all sorts of headaches  	
* R (v. 3.3.0) should be located in the user's PATH  

#### Running analysis  
```git clone https://github.com/SchlossLab/Sze_CRCMetaAnalysis_XXXX_2016.git```  
```make write.paper```
