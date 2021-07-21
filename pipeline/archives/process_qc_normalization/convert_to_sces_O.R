# this script has a function that makes a sce from cell ranger output. 
# NOTE: this is an old script and we don't really use it. Use the one called convert_to_sces.R instead.

# as inputs, it needs a list of filepaths that contains all the data we need. 

library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

convert_to_sces <- function(path_to_10X, id) {

    sce <- read10xCounts(samples = path_to_10X) # depending on which number we are at, pick the data set from the list with the correct index 
    
    # rename reporter genes
    rowData(sce)$Symbol[rowData(sce)$ID == 'ENSG00099999996'] <- 'tdTomato' # DH21, 22, 21_control, 22_control
    rowData(sce)$Symbol[rowData(sce)$ID == 'ENSG00099999997'] <- 'eGFP' # DH21, 22, 21_control, 22_control
    rowData(sce)$Symbol[rowData(sce)$ID == 'ENSG00099999998'] <- 'eGFP' # DH23
    rowData(sce)$Symbol[rowData(sce)$ID == 'ARID1a-Tdtomato'] <- 'tdTomato' # DH21_P5, DH21_P5_control
    rowData(sce)$Symbol[rowData(sce)$ID == 'PIK3CA-GFP'] <- 'eGFP' # DH21_P5, DH21_P5_control
    
    # rename colnames the Barcode ids
    colnames(sce) <- sce$Barcode
    
    # make gene symbol the default rowname (legibility), in cases of collision use both ID and Symbol
    rownames(sce) <- uniquifyFeatureNames(ID = rownames(sce), names = rowData(sce)$Symbol)
    
    # Add a column to colData with the name of the experiment (for groupby later on)
    sce$id <- rep(id, dim(sce)[2])
    
    # check if any of the col names have NA, this causes problems later on. we will change it with 'NA' instead. 
    location <- is.na(names(rowData(sce)))
    names(rowData(sce))[location] <- 'NA'
  
  return(sce) #return the list that contains all the sce objects we made in this function 
  
} # end of function 


