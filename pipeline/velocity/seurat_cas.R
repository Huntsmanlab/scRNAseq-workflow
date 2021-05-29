library(here)
library(loomR)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "make a seurat object with cell type information")

parser$add_argument('--sce_cas_path', metavar='FILE', type='character', help="Path raw sce")
parser$add_argument('--path_to_loom_file', metavar='FILE', type='character', help="output path for loom file")
parser$add_argument('--path_to_barcode_tsv', metavar='FILE', type='character', help="output path for tsv file containing barcode information")

args <- parser$parse_args()

make_seurat_object <- function(sce_cas_path, 
                               path_to_loom_file,
                               path_to_barcode_tsv) {
  
  # read the sce with cell type info 
  sce <- readRDS(sce_cas_path)
  sobject = CreateSeuratObject(counts = counts(sce))
  
  sobject <- FindVariableFeatures(sobject)
  
  # add cell type info 
  sobject$cell_type <- sce$cell_type
  
  # make loom file from the sobject with cell assign info
  seurat_cas.loom <- as.loom(sobject, filename = path_to_loom_file, verbose = FALSE)
  
  # save barcode in tsv
  write.table(as.data.frame(colData(sce_pro_cas)$Barcode), 
              file=path_to_barcode_tsv, quote=FALSE, col.names = NA)
  
} # end of function 

make_seurat_object(sce_cas_path = args$sce_cas_path, 
                   path_to_loom_file = args$path_to_loom_file,
                   path_to_barcode_tsv = args$path_to_barcode_tsv) # call the function 







