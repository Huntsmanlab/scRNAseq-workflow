# this scripts contains the markers that will be used in the create_marker_matrix script to when running cell assign 
library(readxl)

# load marker genes list OUR OWN
our_genes <- read_excel("/huntsman/amunzur/projects/paper/gene_list/lineage_markers.xlsx")
names(our_genes) <- c("non_ciliated_cells", "ciliated_cells", "proliferative_cells", "cluster4", "secretory_cells_like" )

# load marker genes list from NATURE
nature_genes <- read_csv("/huntsman/amunzur/projects/paper/gene_list/natureMed_markers.csv")
names(nature_genes) <- c("unciliated_epithelium", "ciliated_epithelium", "stromal_fibroblasts", "endothelium", "macrophage" )

epithelial_gene_symbols <- c('CLDN3', 'KRT8', 'KRT19', 'WFDC2', 'KLF5', 'SDC4', 'UCA1', 'TACSTD2', 'LINC01541', 'ELF3', 'C1orf186', 'DSP', 'CLDN4', 'PERP', 'KRT18', 'CD9', 'USP53')
lymphocyte_gene_symbols <- c("MS4A1", 'CD79A', 'PTPRC', 'CD19', 'BANK1', 'CD24', 'IGKC', "CD4","CD2","CD3G","CD3D","CD28","CD3E", "CCL5", "STK17B")
endothelial_gene_symbols <- na.omit(nature_genes$endothelium)
macrophage_gene_symbols <- na.omit(nature_genes$macrophage)
stromal_fibroblast_gene_symbols <- na.omit(nature_genes$stromal_fibroblasts)
plasma_gene_symbols <- c('SDC1', 'IGHG1', 'IGHG2', 'CAV1')



