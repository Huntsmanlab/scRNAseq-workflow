# this scripts contains the markers that will be used in the create_marker_matrix script to when running cell assign 

B_gene_symbols <- c('MS4A1', 'CD79A', 'PTPRC', 'CD19', 'BANK1', 'CD24', 'IGKC')
Plasma_gene_symbols <- c('CD79A', 'PTPRC', 'SDC1', 'IGKC', 'IGHG1', 'IGHG2', 'CAV1')
ovarian_cancer <- c('PAX8', 'WT1', 'STMN1', 'MUC16', 'KL6', 'KL7', 'KL8')

# Epithelial_gene_symbols = c("FUT4","EPCAM", "CDH1", "CLDN3", "CLDN4", "KRT8", "KRT18", "KRT19", 
#                             "RPS3", "KRT7", "HMGA1", "TUBB", "TPM4", "S100A16", "MMP7", "MARCKSL1", "MDK", "LCN2", 
#                             'CD24', 'OCLN', 'KRT19', 'DSP', 'CLDN4')

Epithelial_gene_symbols <- c('CD24', "FUT4","EPCAM", "CDH1", "CLDN3", 'COL1A2', 'SCGB1D2', 'TFPI2', 'SERPINA5', 'SLC26A2', 'MGST1', 'SLC39A14', 'HMGCR', 'ENPP3', 'GABRP')

Mesenchymal_gene_symbols <- c('CD44', 'CDH2', 'ITGA5', 'VIM', 'FN1', 'S100A4', 'TNC', 'MMP2', 'ACTA2', 'TWIST1', 'WNT5A', 'SNAI2', 'ZEB1', 'ZEB2')

# Epithelial_ciliated_gene_symbols <- c("FUT4","EPCAM", "CDH1", "CLDN3", "CLDN4", "KRT8", "KRT18", "KRT19", "FOXJ1", "DYDC2", "SNTN", "WDR16", "FAM92B", 
#                                      "C20orf85", "TPPP3", "C9orf24", "OMG", "IGFBP7", "FAM183A", "RSPH1", "C1orf194", "C5orf49", "RP11-356K23.1")

Epithelial_ciliated_gene_symbols <- c('C20orf85', 'TPPP3', 'TUBA1A', 'RSPH1', 'TUBB4B', 'C9orf24', 'IGFBP7', 'C5orf49', 'PIFO', 'OMG', 'FOXJ1', 'DYDC2')

# Unciliated_cells <- c('CXCL14', 'STRA6', 'KLK12', 'HES5', 'RGS16', 'FRZB', 'MUC6', 'ASCL2', 'SCX', 'NUPR1')
Proliferative_cells <- c('MKI67', 'HMGB2', 'CENPF', 'TOP2A', 'CCNB1', 'PCLAF', 'UBE2C', 'TYMS')

# endo_stem_cell <- read.csv('/huntsman/mdouglas/dh_organoid/data/cellassign/endometrialStemCell.csv', header = FALSE)
# mes_stem_cell <- read.csv('/huntsman/mdouglas/dh_organoid/data/cellassign/mesenchymalStemCell.csv', header = FALSE)

# endo_stem_cell <- data.frame(lapply(endo_stem_cell ,as.character), stringsAsFactors=FALSE) %>% unlist(use.names = FALSE)
# mes_stem_cell <- data.frame(lapply(mes_stem_cell ,as.character), stringsAsFactors=FALSE) %>% unlist(use.names = FALSE)
endo_stem_cell <- c('CD34', 'CD44', 'SOX2')
mes_stem_cell <- c('ENG', 'THY1', 'ALCAM')

Cluster2_like <- c('RHBTB3', 'SPTSSB', 'GUCY1A3', 'C6orf15', 'LTF', 'TNFRSF18', 'ECI2', 'OLFM4', 'ALPP', 'SERPINA5', 'COL1A2', 'TFF3', 'MUC5B', 'PTGS2')