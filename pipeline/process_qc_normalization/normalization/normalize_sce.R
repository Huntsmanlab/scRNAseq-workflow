normalize_sce <- function(sce, 
                          output_file_name) {
  
  print("Started normalization.")
  
  # set seed for PCA, TSNE and UMAP
  set.seed(1111)
  
  # load sce
  sce_qc <- sce

  min_size <- min(150, floor(dim(sce_qc)[2] * 0.3))
  max_win <- min(101, min_size + 1)
  
  clusters <- quickCluster(sce_qc, min.size = 50, assay.type="counts", method = "igraph", use.ranks = FALSE, BSPARAM = IrlbaParam())
  sce_qc <- computeSumFactors(sce_qc, assay.type="counts", sizes = seq(21, max_win, 5), min.mean = 0.1, clusters = clusters)
  
  # Ensure that size factors are non-zero and non-negative before normalizing 
  if(summary(sizeFactors(sce_qc))[['Min.']] <= 0){
    stop('Negative size factors, cannot progess')
  }
  
  
  # Normalize using size factors (within batch)
  sce_qc <- logNormCounts(sce_qc)
  
  print("Finished normalization, saving the sce_norm object right now.")
  print(output_file_name)

  # save the output 
  saveRDS(sce_qc, file = output_file_name)
}