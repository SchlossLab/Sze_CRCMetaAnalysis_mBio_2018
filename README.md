## Making Sense of the Noise: Leveraging Existing 16S rRNA Gene Surveys to Identify Key Community Members in Colorectal Tumors


**Background.** An increasing body of literature suggests that both individual and collections of bacteria are associated with the progression of colorectal cancer. As the number of studies investigating these associations increases and the number of subjects in each study increases, a meta-analysis to identify the associations that are the most predictive of disease progression is warranted. For our meta-analysis, we analyzed previously published 16S rRNA gene sequencing data collected from feces (1737 individuals from 8 studies) and colon tissue (492 total samples from 350 individuals from 7 studies).


**Results.** We quantified the odds ratios for individual bacterial genera that were associated with an individual having tumors relative to a normal colon. Among the stool samples, there were no genera that had a significant odds ratio associated with adenoma and there were 8 genera with significant odds ratios associated with carcinoma. Similarly, among the tissue samples, there were no genera that had a significant odds ratio associated with adenoma and there were 3 genera with significant odds ratios associated with carcinoma. Among the significant odds ratios, the association between individual taxa and tumor diagnosis was equal or below 7.11. Because individual taxa had limited association with tumor diagnosis, we trained Random Forest classification models using the genera with the five highest and lowest odds ratios, using the entire collection of genera found in each study, and using operational taxonomic units defined based on a 97% similarity threshold.
All training approaches yielded similar classification success as measured using the Area Under the Curve. The ability to correctly classify individuals with adenomas was poor and the ability to classify individuals with carcinomas was considerably better using sequences from stool or tissue.


**Conclusions.** This meta-analysis confirms previous results indicating that individuals with adenomas cannot be readily classified based on their bacterial community, but that those with carcinomas can. Regardless of the dataset, we found a subset of the fecal community that was associated with carcinomas was as predictive as the full community.





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
