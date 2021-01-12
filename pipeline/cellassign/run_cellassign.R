library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))
source(here('pipeline', 'cellassign', 'create_marker_matrix.R'))
source(here('pipeline', 'cellassign', 'assign_cell_type.R'))

parser <- ArgumentParser(description = "run cell assign on a sample")

parser$add_argument('--sce_clus', metavar='DIRECTORY', type='character', help="path to sce clus")

parser$add_argument('--sce_cas', metavar='DIRECTORY', type='character', help="where we save the marker matrix")

parser$add_argument('--marker_mat_path', metavar='DIRECTORY', type='character', help="where we save the sce with cell type info matrix")

parser$add_argument('--cellassignment_path', metavar='DIRECTORY', type='character', help="where we run cell assignment")

args <- parser$parse_args()

cellAssign <- function(sce_clus, 
                       sce_cas, 
                       marker_mat_path,
                       cellassignment_path){
  
  marker_mat <- create_marker_mat(marker_mat_path) 
                                 
  
  #### Run cell assign
  find_cell_type(path_to_sceqc = sce_clus, 
                 path_to_marker_mat = marker_mat_path,  
                 output_file_name1 = sce_cas, 
                 output_file_name2 = cellassignment_path)
}


cellAssign(sce_clus = args$sce_clus, 
           sce_cas = args$sce_cas, 
           marker_mat_path = args$marker_mat_path,
           cellassignment_path = args$cellassignment_path)


