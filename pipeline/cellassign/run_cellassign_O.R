library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))
source(here('pipeline', 'cellassign', 'create_marker_matrix.R'))
source(here('pipeline', 'cellassign', 'assign_cell_type.R'))

id <- 'VOA11068_ENOC-DH13-DH18'

marker_mat <- create_marker_mat(output_file_name = paste('../data/cellassign', id, 'marker.mat.rds', sep = '/'), 
                                mode = 'split_epithelial') 


#### Run cell assign
find_cell_type(path_to_sceqc = paste('../data/integrated', id, 'integrated_sce.rds', sep = '/'), 
               path_to_marker_mat = paste('../data/cellassign', id, 'marker.mat.rds', sep = '/'),  
               output_file_name1 = paste('../data/cellassign', id, 'sce_norm_cas.rds', sep = '/'), 
               output_file_name2 = paste('../data/cellassign', id, 'cellassignment.rds', sep = '/'))

sce <- SingleCellExperiment(list(counts = GetAssayData(seurat)))
sce$id <- seurat$orig.ident

saveRDS(sce, file = '/huntsman/amunzur/data/integrated/VOA11068_ENOC-DH13-DH18/integrated_sce.rds')
