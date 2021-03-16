suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(cellassign)
  library(org.Hs.eg.db)
  library(biomaRt)
  library(aargh)
})

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

library(readxl)

# load marker genes list OUR OWN
our_genes <- read_excel("/huntsman/amunzur/projects/paper/gene_list/lineage_markers.xlsx")
names(our_genes) <- c("non_ciliated_cells", "ciliated_cells", "proliferative_cells", "cluster4", "secretory_cells_like" )

# load marker genes list from NATURE
nature_genes <- read_csv("/huntsman/amunzur/projects/paper/gene_list/natureMed_markers.csv")
names(nature_genes) <- c("unciliated_epithelium", "ciliated_epithelium", "stromal_fibroblasts", "endothelium", "macrophage" )


create_marker_mat <- function(output_file_name){
  
  # load the marker genes provided by cellassign
  # data(example_TME_markers)
  
  # contruct gene marker matrix
  # marker_list <- example_TME_markers$ensembl
  # marker_list <- lapply(marker_list, function(list) list[! names(list) %in% c('VIM')]) # remove VIM
  
  epithelial_cells <- c('CLDN3', 'KRT8', 'KRT19', 'WFDC2', 'KLF5', 'SDC4', 'UCA1', 'TACSTD2', 'LINC01541', 'ELF3', 'C1orf186', 'DSP', 'CLDN4', 'PERP', 'KRT18', 'CD9', 'USP53')
  lymphocyte_cells <- c("MS4A1", 'CD79A', 'PTPRC', 'CD19', 'BANK1', 'CD24', 'IGKC', "CD4","CD2","CD3G","CD3D","CD28","CD3E", "CCL5", "STK17B")
  endothelial_cells <- na.omit(nature_genes$endothelium)
  macrophage_cells <- na.omit(nature_genes$macrophage)
  stromal_fibroblast_cells <- na.omit(nature_genes$stromal_fibroblasts)
  plasma_cells <- c('SDC1', 'IGHG1', 'IGHG2', 'CAV1')
  
  marker_list <- list(epithelial_cells = c('CLDN3', 'KRT8', 'KRT19', 'WFDC2', 'KLF5', 'SDC4', 'UCA1', 'TACSTD2', 'LINC01541', 'ELF3', 'C1orf186', 'DSP', 'CLDN4', 'PERP', 'KRT18', 'CD9', 'USP53'),
                      lymphocyte_cells = c("MS4A1", 'CD79A', 'PTPRC', 'CD19', 'BANK1', 'CD24', 'IGKC', "CD4","CD2","CD3G","CD3D","CD28","CD3E", "CCL5", "STK17B"),
                      endothelial_cells = na.omit(nature_genes$endothelium),
                      macrophage_cells = na.omit(nature_genes$macrophage),
                      stromal_fibroblast_cells = na.omit(nature_genes$stromal_fibroblasts),
                      plasma_cells = c('SDC1', 'IGHG1', 'IGHG2', 'CAV1'))
  
  gene_symbols = Reduce(union, marker_list)
  
  id_map <- select(org.Hs.eg.db, keys = gene_symbols, columns = c("ENSEMBL"), keytype = "SYMBOL")
  colnames(id_map) <- c("hgnc_symbol", "ensembl_gene_id")
  
  id_map <- id_map[!duplicated(id_map$hgnc_symbol), ] # remove duplicate genes, keep first occurence
  
  cell_types <- c('epithelial_cells', 
                  'lymphocyte_cells', 
                  'endothelial_cells', 
                  'macrophage_cells', 
                  'stromal_fibroblast_cells', 
                  'plasma_cells')
  
  cell_markers <- list(epithelial_cells,
                       lymphocyte_cells,
                       endothelial_cells,
                       macrophage_cells,
                       stromal_fibroblast_cells,
                       plasma_cells)
  
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