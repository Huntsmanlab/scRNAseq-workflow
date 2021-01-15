library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "cluster the sce")

parser$add_argument('--path_to_sce_red', metavar='DIRECTORY', type='character', help="Path normalized sce")

parser$add_argument('--output_file_name', metavar='FILE', type='character', help="Path to normalized sce")

parser$add_argument('--seed', metavar='FILE', type='integer', help="seed for clustering")

parser$add_argument('--k_value', metavar='FILE', type='integer', help="k val for clustering")

args <- parser$parse_args()

i_love_me <- function(path_to_sce_red, 
                      output_file_name, 
                      seed,
                      k_value) {
  
  # load sce
  sce_red <- readRDS(path_to_sce_red)
  
  g <- buildSNNGraph(sce_red, use.dimred = "PCA", k = k_value) # shared neighbor weighting
  cluster <- igraph::cluster_walktrap(g)$membership # extract clustering info from the graph that we just made 
  sce_red$cluster <- factor(cluster) # add that info to the sce 
  
  # save the output 
  saveRDS(sce_red, file = output_file_name)
  
} # end of function

i_love_me(path_to_sce_red = args$path_to_sce_red, 
          output_file_name = args$output_file_name, 
          seed = args$seed, 
          k_value = args$k_value)