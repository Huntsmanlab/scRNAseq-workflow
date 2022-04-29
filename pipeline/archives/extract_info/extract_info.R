library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "extract info from data to prepare for velocity analysis")

parser$add_argument('--path_to_data', metavar='DIRECTORY', type='character', help="Path to sce with clustering info")
parser$add_argument('--counts_matrix', metavar='DIRECTORY', type='character', help="Path to counts matrix")
parser$add_argument('--pca_coords', metavar='DIRECTORY', type='character', help="Path to PCA coordinates")
parser$add_argument('--tsne_coords', metavar='DIRECTORY', type='character', help="Path to t-SNE coordinates")
parser$add_argument('--umap_coords', metavar='DIRECTORY', type='character', help="Path to UMAP coordinates")
parser$add_argument('--clusters', metavar='DIRECTORY', type='character', help="Path to clustering info")
parser$add_argument('--data_type', metavar = 'FILE', type = 'character', help = 'Is this data a sce or seurat?')

args <- parser$parse_args()

extract_velocity_info <- function(path_to_data, 
                                  counts_matrix, 
                                  pca_coords, 
                                  tsne_coords, 
                                  umap_coords, 
                                  clusters,
                                  data_type) {
  
  # load the data, sce or seurat
  object <- readRDS(path_to_data)
  
  if (data_type == "sce"){
    
    # counts matrix
    counts_matrix_data <- as.data.frame(logcounts(object)) %>% 
      rownames_to_column()
    
    # dim reduction
    pca_data <- as.data.frame(object@int_colData@listData[["reducedDims"]]@listData[["PCA"]])
    tsne_data <- as.data.frame(object@int_colData@listData[["reducedDims"]]@listData[["TSNE"]])
    umap_data <- as.data.frame(object@int_colData@listData[["reducedDims"]]@listData[["UMAP"]])
    
    # clusters
    clusters_data <- as.data.frame(object$cluster)
    
  } else {
    
    # integrated counts matrix
    counts_matrix_data <- as.data.frame(GetAssayData(object = object[["integrated"]], slot = "data"))
    
    # dim reduction
    pca_data <- as.data.frame(object@reductions[["pca"]]@cell.embeddings)
    tsne_data <- as.data.frame(object@reductions[["tsne"]]@cell.embeddings)
    umap_data <- as.data.frame(object@reductions[["umap"]]@cell.embeddings)
    
    # clusters
    object$seurat_clusters <- as.numeric(object$seurat_clusters)
    clusters_data <- as.data.frame(object$seurat_clusters)

  }
  
  # save outputted data frames
  write_csv(counts_matrix_data, counts_matrix)
  write_csv(pca_data, pca_coords)
  write_csv(tsne_data, tsne_coords)
  write_csv(umap_data, umap_coords)
  write_csv(clusters_data, clusters)
  
} # end of function

extract_velocity_info(path_to_data = args$path_to_data, 
                      counts_matrix = args$counts_matrix, 
                      pca_coords = args$pca_coords, 
                      tsne_coords = args$tsne_coords, 
                      umap_coords = args$umap_coords, 
                      clusters = args$clusters, 
                      data_type = args$data_type)














