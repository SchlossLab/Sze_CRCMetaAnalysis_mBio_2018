### List of commands to render final pdf
### renders both main text and supplemental

# Load needed functions
source('code/functions.R')

# Load needed library
loadLibs(c("knitr", "rmarkdown"))

# Render the final pdfs
render('submission/manuscript.Rmd', clean = FALSE)

render('submission/supplemental.Rmd', clean=FALSE)
