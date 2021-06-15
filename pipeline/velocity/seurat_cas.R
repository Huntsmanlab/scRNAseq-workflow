library(here)
library(loomR)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "make a seurat object with cell type information")

parser$add_argument('--cell_type_csv', metavar='FILE', type='character', help="Path to cell type csv file")
parser$add_argument('--sce_clus', metavar='FILE', type='character', help='Path to sce clus rds')
parser$add_argument('--path_to_loom_file', metavar='FILE', type='character', help="output path for loom file")

args <- parser$parse_args()

make_seurat_object <- function(cell_type_csv, 
                               sce_clus,
                               path_to_loom_file) {
  
  # add cell type information to sce 
  sce <- readRDS(sce_clus)
  cell_type <- read_csv(cell_type_csv)
  sce$cell_type <- cell_type$cell_type
  
  # create seurat object
  sobject = CreateSeuratObject(counts = counts(sce))
  sobject <- FindVariableFeatures(sobject)
  
  # add cell type info 
  sobject$cell_type <- sce$cell_type
  
  # make loom file from the sobject with cell assign info
  seurat_cas.loom <- as.loom(sobject, filename = path_to_loom_file, verbose = FALSE)
  
  # save barcode in tsv
  #write.table(as.data.frame(colData(sce_pro_cas)$Barcode), 
  #            file=path_to_barcode_tsv, quote=FALSE, col.names = NA)
  
} # end of function 

make_seurat_object(cell_type_csv = args$cell_type_csv, 
                   sce_clus = args$sce_clus,
                   path_to_loom_file = args$path_to_loom_file) 


# cell_type_csv = '../data/cellassign/DH21_NEW/cell_types.csv'
# sce_clus = '../data/clustered/sce/DH21_NEW/sce_clus.rds'
# path_to_loom_file = '../data/loom_cas/DH21_NEW/seurat_cas.loom'
# make_seurat_object(cell_type_csv = cell_type_csv, 
#                    sce_clus = sce_clus,
#                    path_to_loom_file = path_to_loom_file) 