library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "make a seurat object with cell type information")

parser$add_argument('--path_to_sce_raw', metavar='FILE', type='character', help="Path raw sce")
parser$add_argument('--path_to_loom_file', metavar='FILE', type='character', help="output path")

args <- parser$parse_args()

make_seurat_object <- function(path_to_sce_raw, 
                               path_to_loom_file) {
  
  # read the sce with cell type info 
  sce <- readRDS(path_to_sce_raw)
  sobject = CreateSeuratObject(counts = counts(sce))
  
  sobject <- FindVariableFeatures(sobject)
  
  # add cell type info 
  sobject$cell_type <- sce$cell_type
  
  # name the loom file 
  # loom_name = paste(unique(sce$id), 'seurat', sep = '_')
  
  # make loom file from the sobject with cell assign info
  seurat_cas.loom <- as.loom(sobject, filename = path_to_loom_file, verbose = FALSE)
  
} # end of function 

make_seurat_object(path_to_sce_raw = args$path_to_sce_raw, 
                   path_to_loom_file = args$path_to_loom_file) # call the function 







