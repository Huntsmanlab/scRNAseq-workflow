library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))
source(here('pipeline', 'cellassign', 'scripts', 'create_marker_matrix_GC.R'))
source(here('pipeline', 'cellassign', 'scripts', 'assign_cell_type.R'))

sample_id <- 'DH30_GC_wtFOXL2_new'

marker_mat <- create_marker_mat(output_file_name = paste('../data/cellassign', sample_id, 'marker.mat.rds', sep = '/'))
 #                               mode = 'split_epithelial')



#### Run cell assign
find_cell_type(path_to_sce = '/huntsman/harpcheng/data/qc/DH30_GC_wtFOXL2_new/sce_qc.rds', 
               path_to_marker_mat = paste('/huntsman/harpcheng/data/cellassign', sample_id, 'marker.mat.rds', sep = '/'),  
               output_file_name1 = paste('/huntsman/harpcheng/data/cellassign', sample_id, 'sce_cas.rds', sep = '/'), 
               output_file_name2 = paste('/huntsman/harpcheng/data/cellassign', sample_id, 'cellassignment.rds', sep = '/'))

sce <- SingleCellExperiment(list(counts = GetAssayData(seurat)))
sce$id <- seurat$orig.ident

saveRDS(sce, file = '/huntsman/amunzur/data/cellassign/DH17/sce_pro_cas.rds')
