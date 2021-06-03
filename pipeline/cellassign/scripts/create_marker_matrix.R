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


create_marker_mat <- function(sample_type){
  # we build the gene marker matrix based on our sample type: whether the sample is granulosa cells ("GC") or ...

  # marker genes for some common cells:
  epithelial_cells <- c('CLDN3', 'KRT8', 'KRT19', 'WFDC2', 'KLF5', 'SDC4', 'UCA1', 'TACSTD2', 'LINC01541', 'ELF3', 'C1orf186', 'DSP', 'CLDN4', 'PERP', 'KRT18', 'CD9', 'USP53')
  lymphocyte_cells <- c("MS4A1", 'CD79A', 'PTPRC', 'CD19', 'BANK1', 'CD24', 'IGKC', "CD4","CD2","CD3G","CD3D","CD28","CD3E", "CCL5", "STK17B")
  endothelial_cells <- c("ADGRL4", "VWF", "PCDH17", "PECAM1", "RNASE1")
  macrophage_cells <- c("AIF1","MS4A6A","HLA-DQA1", "CD14","IL1B","LYZ","CYBB","HLA-DQB1", "PLEK","HLA-DPA1", "HLA-DPB1", "HLA-DRB5", "HLA-DRB6")
  stromal_fibroblast_cells <- c("FHAD","DCN","COL6A1","CRISPLD2","COL6A3","LUM","COL5A1")
  plasma_cells <- c('SDC1', 'IGHG1', 'IGHG2', 'CAV1')
  granulosa_cells <- c("FSHR", "CYP19A1",	"AMH",	"MGARP", "GLDC", "CHST8", "GPX3",	"FOXL2", "MCAM", "INHA")
  theca_cells <- c("CYP17A1",	"INSL3", "FBLN5",	"OGN", "RAMP2")
  mesothelial_cells <- c("FRAS1",	"RSPO1",	"MSLN",	"LRRN4",	"CALB2",	"KRT5", "MEDAG")
  
  if(sample_type=="GC"){
    marker_list <- list(granulosa_cells=granulosa_cells,
                        theca_cells=theca_cells,
                        mesothelial_cells=mesothelial_cells,
                        stromal_fibroblast_cells=stromal_fibroblast_cells,
                        endothelial_cells=endothelial_cells,
                        plasma_cells=plasma_cells,
                        epithelial_cells=epithelial_cells,
                        lymphocyte_cells=lymphocyte_cells,
                        macrophage_cells=macrophage_cells)
    
  } else { marker_list <- list(epithelial_cells = epithelial_cells,
                              lymphocyte_cells = lymphocyte_cells,
                              endothelial_cells = endothelial_cells,
                              macrophage_cells = macrophage_cells,
                              stromal_fibroblast_cells = stromal_fibroblast_cells,
                              plasma_cells = plasma_cells) }

  
  gene_symbols = Reduce(union, marker_list)
  id_map <- select(org.Hs.eg.db, keys = gene_symbols, columns = c("ENSEMBL"), keytype = "SYMBOL")
  colnames(id_map) <- c("hgnc_symbol", "ensembl_gene_id")
  id_map <- id_map[!duplicated(id_map$hgnc_symbol), ] # remove duplicate genes, keep first occurence
  
  cell_types <- names(marker_list)
  cell_markers <- marker_list[cell_types]
    
  for(i in 1:length(cell_types)){
    map <- id_map[id_map$hgnc_symbol %in% cell_markers[[i]], ]
    marker_list[[cell_types[[i]]]] <- setNames(map$ensembl_gene_id, map$hgnc_symbol)}
  
  # construct the maker gene matrix
  marker_mat <- marker_list_to_mat(marker_list)
  
  return(marker_mat)
}


