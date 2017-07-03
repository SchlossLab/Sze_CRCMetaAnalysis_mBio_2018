### Place to store useful functions that will be used repeatedly throughout
## Marc Sze

#This function loads and installs libraries if they are not already loaded by the user 
loadLibs <- function(deps){
  for (dep in deps){
    if (dep %in% installed.packages()[,"Package"] == FALSE){
      install.packages(as.character(dep), quiet=TRUE);
    }
    library(dep, verbose=FALSE, character.only=TRUE)
  }
}


# This function reads in a triangle distance matrix and converts it into a square
read.dist <- function(file, input='lt', make.square=T, diag=0){
  if(input=='lt'){
    stuff <- scan(file, what='character') #gets all of the elements of the file
    n <- as.numeric(stuff[1]) # n = number of groups/samples in file
    stuff <- stuff[-1] # removes number of groups from list of stuff
    m <- data.frame(matrix(NA, nrow=n, ncol=n) ) #makes empty matrix based on number of groups
    diag(m) <- diag #fills in diagonal with specified value
    
    c <- 1 # c keeps track of postion in stuff vector
    for(i in 1:n){
      group <- stuff[c] #get group name
      colnames(m)[i] <- group
      rownames(m)[i] <- group
      
      if(i > 1){
        m[i,1:(i-1)] <- stuff[(c+1):(c+i-1)] #fills in matrix with values from stuff
      }
      c<-c+i #this because math
    }
    if(make.square){
      m[upper.tri(m)] <- t(m)[upper.tri(m)] #fills in upper triangle
    }
  }
  
  if(input=='square'){ #reads in square matrix
    m<-read.table(file, skip=1, row.names=1)
    colnames(m) <- rownames(m)
  }
  return(m)
}


# Create a label key for the facet wrap function for ggplot2
createTaxaLabeller <- function(taxaTable){
  # collects all entries from the first data list that is not unclassified
  dataList <- apply(taxaTable, 1, function(x) x[x != "unclassified"])
  # creates a vector of the lowest ID taxonomy
  if (class(dataList) == "list"){
    
    tempCall <- unname(unlist(lapply(dataList, function(x) x[length(x)])))
    
  } else{
    
    tempCall <- apply(dataList, 2, function(x) x[length(x[x != "Bacteria"])+1])
  }
  
  # assigns names to the vector that are the OTU labels
  names(tempCall) <- rownames(taxaTable)
  # returns that vector
  return(tempCall)
}

# Creates names for the get_tax_level_shared function
get_tax_substring <- function(tax, tax_level){
  substring <- unlist(strsplit(tax, ";"))[tax_level]
  paste(substring, collapse='.')
}
get_tax_name <- function(tax_file, tax_level){
  
  
  tax_data <- read.table(file=tax_file, header=T, stringsAsFactors=F)
  taxonomy <- tax_data$Taxonomy
  taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
  taxonomy <- gsub('"', '', taxonomy)
  
  tax_substring <- sapply(taxonomy, get_tax_substring, tax_level)
  
  names(tax_substring) <- tax_data$OTU
  
  tax_substring
}

# Get the total number based on shared file and tax file.
get_tax_level_shared <- function(shared_file, tax_file, tax_level){
  
  shared_otus <- shared_file
  is_present <- apply(shared_otus, 2, sum) > 0
  shared <- shared_otus[,is_present]
  
  taxonomy <- get_tax_name(tax_file, tax_level)
  taxonomy <- taxonomy[colnames(shared)]
  unique_taxa <- levels(as.factor(taxonomy))
  
  shared_tax_level <- NULL
  
  for(ut in unique_taxa){
    otus <- names(taxonomy[taxonomy %in% ut])
    sub_shared <- shared_otus[,colnames(shared_otus) %in% otus]
    
    if(is.null(dim(sub_shared))){
      shared_tax_level <- cbind(shared_tax_level, sub_shared)
    } else {
      tax_level_count <- apply(sub_shared, 1, sum)
      shared_tax_level <- cbind(shared_tax_level, tax_level_count)
    }
  }
  colnames(shared_tax_level) <- unique_taxa
  rownames(shared_tax_level) <- rownames(shared)
  return(shared_tax_level)
}




















