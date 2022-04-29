suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(cellassign)
  library(org.Hs.eg.db)
  library(biomaRt)
  library(aargh)
})

source(here('pipeline', 'cellassign', 'marker_genes_list.R'))

######################################
#### Build the CellAssign Marker Matrix
######################################
# Inputs:
# 
# Outputs:
#
#
# To see/save the list of genes that are used as markers in a human-readable way:
# marks <- read.table('/huntsman/mdouglas/dh_organoid/data/cellassign/allCellMarkers.csv', 
#                     sep = ",", header = FALSE, col.names = paste0("V",seq_len(25)), fill = TRUE)
# marks <- t(data.frame(lapply(marks,as.character), stringsAsFactors=FALSE))
# png("allCellMarkers.png", height = 50*nrow(marks), width = 450*ncol(summary_DA))
# grid.table(marks)
# dev.off()

create_marker_mat <- function(output_file_name, 
                              mode){
  
  # load the marker genes provided by cellassign
  data(example_TME_markers)
  
  # contruct gene marker matrix
  marker_list <- example_TME_markers$ensembl
  marker_list <- lapply(marker_list, function(list) list[! names(list) %in% c('VIM')]) # remove VIM
  
  if(mode == "split_epithelial"){
    
    gene_symbols = Reduce(union, list(B_gene_symbols, 
                                      Plasma_gene_symbols, 
                                      Epithelial_gene_symbols, 
                                      Epithelial_ciliated_gene_symbols, 
                                      endo_stem_cell, 
                                      mes_stem_cell, 
                                      ovarian_cancer, 
                                      Mesenchymal_gene_symbols, 
                                      Proliferative_cells))
    
    id_map <- select(org.Hs.eg.db, keys = gene_symbols, columns = c("ENSEMBL"), keytype = "SYMBOL")
    colnames(id_map) <- c("hgnc_symbol", "ensembl_gene_id")

    id_map <- id_map[!duplicated(id_map$hgnc_symbol), ] # remove duplicate genes, keep first occurence
    
    cell_types <- c('B cells', 
                    'Plasma cells', 
                    'Epithelial cells', 
                    'Epithelial ciliated cells', 
                    'Endometrial stem cells', 
                    'Mesenchymal stem cells', 
                    'Mesenchymal cells', 
                    'Proliferative cells')
    
    cell_markers <- list(B_gene_symbols, 
                         Plasma_gene_symbols, 
                         Epithelial_gene_symbols, 
                         Epithelial_ciliated_gene_symbols, 
                         endo_stem_cell, 
                         mes_stem_cell, 
                         ovarian_cancer, 
                         Mesenchymal_gene_symbols, 
                         Proliferative_cells)
    
    
  } else if(mode == 'nosplit_epithelial'){
    
    gene_symbols = Reduce(union, list(B_gene_symbols, 
                                      Plasma_gene_symbols, 
                                      Epithelial_gene_symbols, 
                                      endo_stem_cell, 
                                      mes_stem_cell, 
                                      ovarian_cancer, 
                                      Mesenchymal_gene_symbols, 
                                      Proliferative_cells))
    
    id_map <- select(org.Hs.eg.db, keys = gene_symbols, columns = c("ENSEMBL"), keytype = "SYMBOL")
    colnames(id_map) <- c("hgnc_symbol", "ensembl_gene_id")
    
    id_map <- id_map[!duplicated(id_map$hgnc_symbol), ] # remove duplicate genes, keep first occurence
    
    cell_types <- c('B cells', 
                    'Plasma cells', 
                    'Epithelial cells', 
                    'Endometrial stem cells', 
                    'Mesenchymal stem cells', 
                    'Mesenchymal cells', 
                    'Proliferative cells')
    
    cell_markers <- list(B_gene_symbols, 
                         Plasma_gene_symbols, 
                         Epithelial_gene_symbols, 
                         endo_stem_cell, 
                         mes_stem_cell, 
                         ovarian_cancer, 
                         Mesenchymal_gene_symbols, 
                         Proliferative_cells)
  }
  
  modify_marker_list <- function(cell_types, cell_markers, marker_list, id_map){
    
    for(i in seq(1, length(cell_types), 1)){
      map <- id_map[id_map$hgnc_symbol %in% cell_markers[[i]], ]
      marker_list[[cell_types[[i]]]] <- setNames(map$ensembl_gene_id, map$hgnc_symbol)
    }
    
    return(marker_list)
  }
  
  marker_list <- modify_marker_list(cell_types, cell_markers, marker_list, id_map)
  
  # construct the maker gene matrix
  marker_mat <- marker_list_to_mat(marker_list)
  
  saveRDS(marker_mat, file = output_file_name)
  
  return(marker_mat)
  
}
