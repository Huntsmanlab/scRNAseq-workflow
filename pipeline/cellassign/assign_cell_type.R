suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(tidyverse)
  library(aargh)
  library(cellassign)
})

find_cell_type <- function(path_to_sceqc, 
                           path_to_marker_mat,  
                           output_file_name1, 
                           output_file_name2){
  
  # size factors already computed in qc step
  sce_qc <- readRDS(path_to_sceqc)
  old_rownames <- rownames(sce_qc)
  rownames(sce_qc) <- rowData(sce_qc)$ID
  
  # construct the maker gene matrix
  marker_mat <- readRDS(path_to_marker_mat)
  
  sce_marker <- sce_qc[intersect(rownames(marker_mat), rownames(sce_qc)), ]
  marker_mat <- marker_mat[intersect(rownames(marker_mat), rownames(sce_qc)), ] 
  
  cas <- cellassign(exprs_obj = sce_marker,
                    marker_gene_info = marker_mat,
                    s = sizeFactors(sce_marker))
  
  sce_qc$cell_type <- cas$cell_type
  rownames(sce_qc) <- old_rownames
  
  saveRDS(sce_qc, file = output_file_name1)
  saveRDS(cas, file = output_file_name2)
  
}
