library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))
source(here('pipeline', 'cellassign', 'scripts', 'create_marker_matrix.R'))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(tidyverse)
  library(aargh)
  library(cellassign)
})

# RUN cell type assignment: save sce with cell type information, save csv with cell type information;
# function returns sce which is used in cell assign report.

run_cellassign <- function(sce_clus, sample_type, cell_type_csv){

  # load marker gene matrix; specify the type of your sample: if granulosa cells, supply "GC"; if..., ...
  marker_mat <- create_marker_mat(sample_type) 
  
  # size factors already computed in qc step
  sce <- readRDS(sce_clus)
  old_rownames <- rownames(sce)
  rownames(sce) <- rowData(sce)$ID
  
  # construct the maker gene matrix
  sce_marker <- sce[intersect(rownames(marker_mat), rownames(sce)), ]
  marker_mat <- marker_mat[intersect(rownames(marker_mat), rownames(sce)), ] 
  
  # assign cell type
  cas <- cellassign(exprs_obj = sce_marker,
                    marker_gene_info = marker_mat,
                    s = sizeFactors(sce_marker))
  
  # make a df with the celltype information 
  cell_type <- as.vector(cas$cell_type) 
  barcodes <- as.vector(colnames(sce)) 
  df <- as.data.frame(cbind(cell_type, barcodes))
  names(df) <- c("cell_type", "barcodes")
  
  # add celltype information to sce
  sce$cell_type <- cas$cell_type
  rownames(sce) <- old_rownames
  
  # save sce with cell type information
  # saveRDS(sce, file = sce_cas_path) 
  
  # save celltype information in a csv file
  dir.create(dirname(cell_type_csv))
  write_csv(df, cell_type_csv) 
  
  # retuen sce with cell type information 
  return(sce)
}

  





