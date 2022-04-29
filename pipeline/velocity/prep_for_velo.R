library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "compute sum factors on a raw sce")

parser$add_argument('--path_to_sce_raw', metavar='FILE', type='character', help="Path raw sce")
parser$add_argument('--path_to_sce_raw_with_sum_factors', metavar='FILE', type='character', help="output path")

args <- parser$parse_args()

prep_for_velocity <- function(path_to_sce_raw, # input 
                              path_to_sce_raw_with_sum_factors){ # output 
  
  # load the processed data - must be before qc 
  sce <- readRDS(path_to_sce_raw) # load the uncorrected object 
  sce <- computeSumFactors(sce)
  saveRDS(sce, path_to_sce_raw_with_sum_factors) # load the uncorrected object with sum factors 

}

prep_for_velocity(path_to_sce_raw = args$path_to_sce_raw, 
                  path_to_sce_raw_with_sum_factors = args$path_to_sce_raw_with_sum_factors)
            

