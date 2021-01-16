# source some packages and general functions
print("Sourcing files.")

library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

# source the files that have the functions for dim red and clustering steps 
source(here('pipeline', 'dimred_cluster', 'dim_reduction', 'dim_reduction.R')) # dim reduction
source(here('pipeline', 'dimred_cluster', 'clustering', 'separate', 'clustering_main.R')) # clustering

parser <- ArgumentParser(description = "Dimensionality reduction and clustering")

# DIM REDUCTION 
parser$add_argument('--path_to_sce_norm', metavar='DIRECTORY', type='character', help="Path normalized sce")
parser$add_argument('--HVG', metavar = 'FILE', type = 'character', help = 'Do you want to subset to HVG?')
parser$add_argument('--top_HVG', metavar='FILE', type='integer', help="percentage for HVG")
parser$add_argument('--top_PCs', metavar='FILE', type='character', help="number PCs to used")

# CLUSTERING
parser$add_argument('--output_file_name', metavar='FILE', type='character', help="Path to normalized sce")
parser$add_argument('--k_value', metavar='FILE', type='integer', help="k val for clustering")

args <- parser$parse_args()

# read the sce, and check if dim red info is already there. if so, skip dim reduction and do clustering only. 
sce <- readRDS(args$path_to_sce_norm)

# even if one of the dim red is missing, repeat dim red step 
if (is.null(sce@int_colData@listData[["reducedDims"]]@listData[["PCA"]]) || # this means or 
    is.null(sce@int_colData@listData[["reducedDims"]]@listData[["TSNE"]]) ||
    is.null(sce@int_colData@listData[["reducedDims"]]@listData[["UMAP"]])) {
  
  print("Missing dim red info. Dim reduction step started.")
  
  sce <- perform_dim_reduction(path_to_sce_norm = args$path_to_sce_norm, 
                               HVG = args$HVG, 
                               top_HVG = args$top_HVG, 
                               top_PCs = args$top_PCs)
  
  cluster_cells(sce = sce,
                output_file_name = args$output_file_name, 
                k_value = args$k_value)
  
  # if all dimred are TRUE (they exist)
} else if (!is.null(sce@int_colData@listData[["reducedDims"]]@listData[["PCA"]]) &&
           !is.null(sce@int_colData@listData[["reducedDims"]]@listData[["TSNE"]]) &&
           !is.null(sce@int_colData@listData[["reducedDims"]]@listData[["UMAP"]])) {
  
  print("It seems the dim red information has been computed before. Only clustering algorithm will be run.")
  
  # skip dim red, do clustering instead
  cluster_cells(sce = sce, # we already loaded sce here 
                output_file_name = args$output_file_name, 
                k_value = args$k_value)
  
}