cluster_cells <- function(sce, 
                          output_file_name, 
                          k_value) {
  
  # load sce
  sce_norm <- sce
  
  print("Started clustering the sce.")  
  g <- buildSNNGraph(sce_norm, use.dimred = "PCA", k = k_value) # shared neighbor weighting
  cluster <- igraph::cluster_walktrap(g)$membership # extract clustering info from the graph that we just made 
  sce_norm$cluster <- factor(cluster) # add that info to the sce 
  
  print("Finished clustering sce.")
  
  # save the output 
  print("Started saving the sce object.")
  saveRDS(sce_norm, file = output_file_name)
  
}

