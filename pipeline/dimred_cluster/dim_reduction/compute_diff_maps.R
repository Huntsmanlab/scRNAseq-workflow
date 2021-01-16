library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "compute diffusion map coordinates")

parser$add_argument('--sce_clus', metavar='DIRECTORY', type='character', help="path to sce clus")
parser$add_argument('--dm_path', metavar='DIRECTORY', type='character', help="where we will save the diffusion map coordinates")

args <- parser$parse_args()

compute_diff_maps <- function(sce_clus, dm_path){
  
  message("Computing diffusion map coordinates. Why not have a cup of tea while waiting?")
  
  sce_clus <- readRDS(sce_clus)
  
  sce_logcounts <- logcounts(sce_clus)  # access log-transformed counts matrix
  sce_logcounts <- as.matrix(sce_logcounts)

  dm <- DiffusionMap(t(sce_logcounts), n_pcs = 30) # outputs an S4 object 
  
  saveRDS(dm, dm_path)
  
  # add diff map information to sce
  tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                    DC2 = eigenvectors(dm)[, 2],
                    DC3 = eigenvectors(dm)[, 3])
  
  sce@int_colData@listData[["reducedDims"]]@listData[["DM"]] <- tmp
  
  

  
  }

compute_diff_maps(sce_clus = args$sce_clus, 
                  dm_path = args$dm_path)


