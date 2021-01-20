library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))
source(here('pipeline', 'cellassign', 'scripts', 'create_marker_matrix.R'))
source(here('pipeline', 'cellassign', 'scripts', 'assign_cell_type.R'))

parser <- ArgumentParser(description = "run cell assign on a sample")

parser$add_argument('--sce_clus', metavar='DIRECTORY', type='character', help="path to sce clus")

parser$add_argument('--cell_type_csv', metavar='DIRECTORY', type='character', help="where we save the cell type information as a csv file")

parser$add_argument('--marker_mat_path', metavar='DIRECTORY', type='character', help="where we save the sce with cell type info matrix")

parser$add_argument('--cellassignment_path', metavar='DIRECTORY', type='character', help="where we run cell assignment")

args <- parser$parse_args()

cellAssign <- function(sce_clus, 
                       marker_mat_path,
                       cell_type_csv,
                       cellassignment_path){
  
  marker_mat <- create_marker_mat(marker_mat_path) 
                                 
  
  #### Run cell assign
  find_cell_type(path_to_sce = sce_clus, 
                 path_to_marker_mat = marker_mat_path,  
                 output_file_name1 = cell_type_csv, 
                 output_file_name2 = cellassignment_path)
}


cellAssign(sce_clus = args$sce_clus, 
           marker_mat_path = args$marker_mat_path,
           cell_type_csv = args$cell_type_csv,
           cellassignment_path = args$cellassignment_path)




