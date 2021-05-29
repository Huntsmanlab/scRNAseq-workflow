convert_to_sces <- function(path_to_10X, id) {
  
  print("Started converting 10X output to an sce object.")
  
  # sometimes our raw data is in a directory called 'filtered_feature_bc_matrix' and sometimes it isn't. 
  # here we account for that: 
  # first check if this directory exists: 
  if(dir.exists(paste(path_to_10X, 'filtered_feature_bc_matrix', sep = '/')) == TRUE) {
    
    path_to_10X <- paste(path_to_10X, 'filtered_feature_bc_matrix', sep = '/')
    
  }
  
  sce <- read10xCounts(samples = path_to_10X) 
  
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
  
  # calculate the QC metrics here. this is QC metrics BEFORE eliminating any low quality cells. 
  # we do elimination in the qc script. 
  # regardless of the method we choose to remove the low quality cells, we need to do this step. 
  # addPerCellQC calculates the QC metrics and adds the data to the colData of our sce. 
  
  if ( grepl("^ENSMUS", sce@rowRanges@elementMetadata@listData[["ID"]][1]) ){               # For mouse
    mt_genes <- grepl("^mt", rowData(sce)$Symbol)
    ribo_genes <- grepl("^Rp[sl][[:digit:]]", rowData(sce)$Symbol)
  } else {                                                                                  # For human
    mt_genes <- grepl("^MT-", rowData(sce)$Symbol)
    ribo_genes <- grepl("^RP[LS]", rowData(sce)$Symbol)
  }
  
  feature_ctrls <- list(mito = rownames(sce)[mt_genes],
                        ribo = rownames(sce)[ribo_genes])
  
  sce <- addPerCellQC(sce, subsets = feature_ctrls)
  
  print("Finished converting 10X output to an sce object.")
  
  return(sce) # return the sce object to be used in the next step
  
} # end of function
