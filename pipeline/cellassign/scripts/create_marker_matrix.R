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


create_marker_mat <- function(cell_types){
  # we build the gene marker matrix based on our sample type: whether the sample is granulosa cells ("GC") or endometrial tissue

  # marker genes for some common cell types (i think this list is derived from https://www.nature.com/articles/s41591-020-1040-z#Sec1 - for endometrium tissue)
  epithelial_cells <- c('CLDN3', 'KRT8', 'KRT19', 'WFDC2', 'KLF5', 'SDC4', 'UCA1', 'TACSTD2', 'LINC01541', 'ELF3', 'C1orf186', 'DSP', 'CLDN4', 'PERP', 'KRT18', 'CD9', 'USP53')
  lymphocyte_cells <- c("MS4A1", 'CD79A', 'PTPRC', 'CD19', 'BANK1', 'CD24', 'IGKC', "CD4","CD2","CD3G","CD3D","CD28","CD3E", "CCL5", "STK17B")
  endothelial_cells <- c("ADGRL4", "VWF", "PCDH17", "PECAM1", "RNASE1")
  macrophage_cells <- c("AIF1","MS4A6A","HLA-DQA1", "CD14","IL1B","LYZ","CYBB","HLA-DQB1", "PLEK","HLA-DPA1", "HLA-DPB1", "HLA-DRB5", "HLA-DRB6")
  stromal_fibroblast_cells <- c("DCN","COL6A1","CRISPLD2","COL6A3","LUM","COL5A1")
  plasma_cells <- c('SDC1', 'IGHG1', 'IGHG2', 'CAV1')

  # for granulosa cells:
  granulosa_cells <- c("FSHR", "CYP19A1",	"AMH",	"MGARP", "GLDC", "CHST8", "GPX3",	"FOXL2", "MCAM", "INHA")
  theca_cells <- c("CYP17A1",	"INSL3", "FBLN5",	"OGN", "RAMP2")
  mesothelial_cells <- c("FRAS1",	"RSPO1",	"MSLN",	"LRRN4",	"CALB2")
  
  # from allCellMarkers.csv:
  Blymphocyte <- c("MS4A1", "CD79A", "PTPRC", "CD19", "BANK1", "CD24", "IGKC")
  Tlymphocyte <- c("CD4","PTPRC","CD2","CD3G","CD3D","CD28","CD3E")
  CytotoxicT <- c("PTPRC","NKG7","CD2","GZMA","CD8A","CD3G","CD3D","CD28","PRF1","CD3E")
  Monocyte <- c("CD4","PTPRC","LYZ","CD33","ITGAX","ITGAM","CD14","FCGR3A")
  Epithelial <- c("FUT4","EPCAM","CDH1","CLDN3","CLDN4","KRT8","KRT18","KRT19","RPS3","KRT7","HMGA1","TUBB","TPM4","S100A16","MMP7","MARCKSL1","MDK","LCN2")
  Myofibroblast <- c("ACTA2","COL1A1","SERPINH1","COL3A1")
  VascularSmoothMuscle<- c("MYLK","MCAM","ACTA2","COL1A1","MYH11","SERPINH1","COL3A1","PLN")
  Endothelial <- c("MCAM","VWF","SERPINH1","EMCN","CLEC14A","CDH5","PECAM1")
  PlasmaB <- c("CD79A", "PTPRC", "SDC1", "IGKC", "IGHG1", "IGHG2", "CAV1")
  EpithelialCiliated <- c("FUT4","EPCAM","CDH1","CLDN3","CLDN4","KRT8","KRT18","KRT19","FOXJ1","DYDC2","SNTN","WDR16","FAM92B","C20orf85","TPPP3","C9orf24","OMG","IGFBP7","FAM183A","RSPH1","C1orf194","C5orf49","RP11-356K23.1")
  EndometrialStem <- c("ABCG2","ENG","VCAM1","KIT","PROM1","CD140b","MCAM","CD166","ITGB1","PECAM1","CD34","CD44","PTPRC","ITGA4","ICAM1","NT5E","THY1","IPO13","LGR5","MSI1","NES","POU5F1","SOX2","SSEA4")
  MesenchymalStem <- c("ENG","THY1","ALCAM","CD34","PTPRC","NT5E")
  Ovarian_cancer <- c("PAX8", "WT1", "STMN1", "MUC16", "KL6", "KL7", "KL8")
  Mesenchymal_cells <- c("CD44", "CDH2", "ITGA5", "VIM", "FN1", "S100A4", "TNC", "MMP2", "ACTA2", "TWIST1", "WNT5A", "SNAI2", "ZEB1", "ZEB2")
  
  # new list:
  epithelial <- c("EPCAM", "CDH1","CLDN3","CLDN4","KRT7","KRT8","KRT18","KRT19","FUT4")
  epithelial_secretory <- c("EPCAM", "KRT7","KRT8","KRT18","PAX8","OVGP1", "KRT19", "LCN2",'WFDC2', 'KLF5', 'SDC4', 'UCA1', 'TACSTD2', 'LINC01541', 'ELF3', 'C1orf186', 'DSP', 'CLDN4', 'PERP', 'KRT18', 'CD9', 'USP53')
  epithelial_ciliated <- c("EPCAM", "KRT7","KRT8","KRT18","FOXJ1", "CAPS", "DYDC2", "CTH", "TP73", "SNTN","WDR16","FAM92B","C20orf85","C20orf88","TPPP3","C9orf24","FAM183A","C1orf194","PIFO")
  endothelial <- c('VWF','PECAM1','MCAM','EMCN','CLEC14A','CDH5')
  stromal_fibroblast <- c("DCN","COL6A1","CRISPLD2","COL6A3","LUM","COL5A1")
  vasular_smooth_muscle <- c("MYLK","ACTA2","MYH11","PLN","SMTN","CNN1")
  lymphocyte <- c("MS4A1", 'CD79A', 'PTPRC', 'CD19', 'BANK1', 'CD24', 'IGKC', "CD4","CD2","CD3G","CD3D","CD28","CD3E", "CCL5", "STK17B")
  mesenchymal <- c("CD44", "CDH2", "ITGA5", "VIM", "FN1", "S100A4", "TNC", "MMP2", "ACTA2", "TWIST1", "WNT5A", "SNAI2", "ZEB1", "ZEB2")
  mesenchymal_stem <- c("ENG","THY1","ALCAM","NT5E")
  endometrial_stem <- c("ABCG2","ENG","VCAM1","KIT","PROM1","CD140b","MCAM","ALCAM","ITGB1","PECAM1","CD34","CD44","PTPRC","ITGA4","ICAM1","NT5E","THY1","IPO13","LGR5","MSI1","NES","POU5F1","SOX2","SSEA4")
  ovarian_cancer <- c("PAX8", "WT1", "MUC16", "KL6", "KL7", "KL8")
  
  # a list of all markers
  marker_list_all <- list(epithelial_cells=epithelial_cells, lymphocyte_cells=lymphocyte_cells, endothelial_cells=endothelial_cells,
                          macrophage_cells=macrophage_cells, stromal_fibroblast_cells=stromal_fibroblast_cells, plasma_cells=plasma_cells,
                          granulosa_cells=granulosa_cells, theca_cells=theca_cells, mesothelial_cells=mesothelial_cells,
                          Blymphocyte=Blymphocyte, Tlymphocyte=Tlymphocyte, CytotoxicT=CytotoxicT, Monocyte=Monocyte,
                          Epithelial=Epithelial, Myofibroblast=Myofibroblast, VascularSmoothMuscle=VascularSmoothMuscle, 
                          Endothelial=Endothelial, PlasmaB=PlasmaB, EpithelialCiliated=EpithelialCiliated, EndometrialStem=EndometrialStem,
                          MesenchymalStem=MesenchymalStem, Ovarian_cancer=Ovarian_cancer, Mesenchymal_cells=Mesenchymal_cells,
                          epithelial=epithelial, epithelial_secretory=epithelial_secretory, epithelial_ciliated=epithelial_ciliated, 
                          endothelial=endothelial, stromal_fibroblast=stromal_fibroblast, vasular_smooth_muscle=vasular_smooth_muscle, 
                          lymphocyte=lymphocyte, mesenchymal=mesenchymal, mesenchymal_stem=mesenchymal_stem, endometrial_stem=endometrial_stem,
                          ovarian_cancer=ovarian_cancer)
  
  # subset to markers you want to include 
  marker_list <- marker_list_all[cell_types]
  
  # if(sample_type=="GC"){ # for granulosa cells
  #   marker_list <- list(granulosa_cells=granulosa_cells,
  #                       theca_cells=theca_cells,
  #                       mesothelial_cells=mesothelial_cells,
  #                       stromal_fibroblast_cells=stromal_fibroblast_cells,
  #                       endothelial_cells=endothelial_cells,
  #                       plasma_cells=plasma_cells,
  #                       lymphocyte_cells=lymphocyte_cells,
  #                       macrophage_cells=macrophage_cells)
  #   
  # } else if (sample_type=="old") { # for endometrial tissue
  #   marker_list <- list(epi_cells = epi_cells,
  #                       lymphocyte_cells = lymphocyte_cells,
  #                       endothelial_cells = endothelial_cells,
  #                       macrophage_cells = macrophage_cells,
  #                       stromal_fibroblast_cells = stromal_fibroblast_cells,
  #                       plasma_cells = plasma_cells,
  #                       epithelial_ciliated = epithelial_ciliated) 
  #   
  # } else if (sample_type=="new") {
  #   marker_list <- list(epithelial_secretory=epithelial_secretory,
  #                       epithelial_ciliated=epithelial_ciliated,
  #                       endothelial=endothelial,
  #                       stromal_fibroblast=stromal_fibroblast,
  #                       vasular_smooth_muscle=vasular_smooth_muscle,
  #                       lymphocyte=lymphocyte,
  #                       #mesenchymal=mesenchymal,
  #                       mesenchymal_stem=mesenchymal_stem,
  #                       endometrial_stem=endometrial_stem)
  #                       #ovarian_cancer=ovarian_cancer)
  #   
  # } else { # from allCellMarkers.csv (also endometrial tissue)
  #     marker_list <- list(Blymphocyte = Blymphocyte,
  #                         Tlymphocyte = Tlymphocyte,
  #                         CytotoxicT = CytotoxicT,
  #                         Monocyte = Monocyte,
  #                         Epithelial = Epithelial,
  #                         Myofibroblast = Myofibroblast,
  #                         VascularSmoothMuscle = VascularSmoothMuscle,
  #                         Endothelial = Endothelial,
  #                         PlasmaB = PlasmaB,
  #                         EpithelialCiliated = EpithelialCiliated,
  #                         EndometrialStem = EndometrialStem,
  #                         MesenchymalStem = MesenchymalStem,
  #                         Ovarian_cancer = Ovarian_cancer,
  #                         Mesenchymal_cells = Mesenchymal_cells)
  #   }
  
  gene_symbols = Reduce(union, marker_list)
  id_map <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, columns = c("ENSEMBL"), keytype = "SYMBOL")
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









