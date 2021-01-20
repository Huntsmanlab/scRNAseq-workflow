library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))
source(here('pipeline', 'cellassign', 'create_marker_matrix.R'))
source(here('pipeline', 'cellassign', 'assign_cell_type.R'))

sample_id <- 'DH15'

marker_mat <- create_marker_mat(output_file_name = paste('../data/cellassign', sample_id, 'marker.mat.rds', sep = '/'), 
                                mode = 'split_epithelial')



#### Run cell assign
find_cell_type(path_to_sceqc = '/huntsman/amunzur/data/processed/DH15/sce_with_sum_factors.rds', 
               path_to_marker_mat = paste('../data/cellassign', sample_id, 'marker.mat.rds', sep = '/'),  
               output_file_name1 = paste('../data/cellassign', sample_id, 'sce_pro_cas.rds', sep = '/'), 
               output_file_name2 = paste('../data/cellassign', sample_id, 'cellassignment.rds', sep = '/'))

sce <- SingleCellExperiment(list(counts = GetAssayData(seurat)))
sce$id <- seurat$orig.ident

saveRDS(sce, file = '/huntsman/amunzur/data/cellassign/DH17/sce_pro_cas.rds')
