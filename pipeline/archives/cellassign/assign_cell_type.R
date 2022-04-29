suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(tidyverse)
  library(aargh)
  library(cellassign)
})

find_cell_type <- function(path_to_sce, 
                           path_to_marker_mat,  
                           output_file_name1, 
                           output_file_name3){

  # size factors already computed in qc step
  sce <- readRDS(path_to_sce)
  old_rownames <- rownames(sce)
  rownames(sce) <- rowData(sce)$ID
  
  # construct the maker gene matrix
  marker_mat <- readRDS(path_to_marker_mat)
  
  sce_marker <- sce[intersect(rownames(marker_mat), rownames(sce)), ]
  marker_mat <- marker_mat[intersect(rownames(marker_mat), rownames(sce)), ] 

  cas <- cellassign(exprs_obj = sce_marker,
                    marker_gene_info = marker_mat,
                    s = sizeFactors(sce_marker))
  
  # make a df with the celltype information 
  cell_type <- as.vector(cas$cell_type) # a vetor with cell type information
  barcodes <- as.vector(colnames(sce)) # a vector with barcodes
  df <- as.data.frame(cbind(cell_type, barcodes))
  names(df) <- c("cell_type", "barcodes") # name the df 
  
  sce$cell_type <- cas$cell_type
  rownames(sce) <- old_rownames
  
  saveRDS(sce, file = output_file_name1) # save sce with cell type information
  write_csv(df, output_file_name3) # save celltype information in a csv file

}



