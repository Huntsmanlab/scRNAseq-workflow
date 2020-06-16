library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "normalize sce after QC")

parser$add_argument('--path_to_QCed_sce', metavar='DIRECTORY', type='character', help="Path QCed sce")

parser$add_argument('--output_file_name', metavar='FILE', type='character', help="Path to normalized sce")

parser$add_argument('--seed', metavar='FILE', type='integer', help="seed for UMAP, TNSE, PCA approximation")

parser$add_argument('--HVG', metavar = 'FILE', type = 'character', help = 'Do you want to subset to HVG?')

parser$add_argument('--top_HVG', metavar='FILE', type='integer', help="percentage for HVG")

args <- parser$parse_args()

normalize_sce <- function(path_to_QCed_sce, 
                          output_file_name, 
                          seed, 
                          HVG, 
                          top_HVG) {
  
  # set seed for PCA, TSNE and UMAP
  set.seed(seed)
  
  # load sce
  sce_qc <- readRDS(path_to_QCed_sce)

  min_size <- min(150, floor(dim(sce_qc)[2] * 0.3))
  max_win <- min(101, min_size + 1)
  
  clusters <- quickCluster(sce_qc, assay.type="counts", method = "igraph", use.ranks = FALSE, BSPARAM = IrlbaParam())
  sce_qc <- computeSumFactors(sce_qc, assay.type="counts", sizes = seq(21, max_win, 5), min.mean = 0.1, clusters = clusters)
  
  # Ensure that size factors are non-zero and non-negative before normalizing 
  if(summary(sizeFactors(sce_qc))[['Min.']] <= 0){
    stop('Negative size factors, cannot progess')
  }
  
  # Normalize using size factors (within batch)
  sce_qc <- logNormCounts(sce_qc)
  
  # deal with HVG
  if (HVG == "yes"){
    
    # pick the HVG
    gene_var <- modelGeneVar(sce_qc)
    
    fraction_HGV <- top_HVG / 100
    chosen <- getTopHVGs(gene_var, prop = top_HVG)

    # subset to HVG
    sce_qc_hvg <- sce_qc[chosen,]

    # add the original count matrix with non HVG to our new sce
    altExp(sce_qc_hvg, "original") <- sce_qc
    altExpNames(sce_qc_hvg)
    
    # to recover original counts 
    # sce_qc_original <- altExp(sce_qc_hvg, "original", withColData=TRUE)
    
  } # end of HVG
  

  # Run UMAP, TSNE, PCA for visualization downstream steps: clustering, cellassign, ect. and save the data to sce 
  sce_qc <- runPCA(sce_qc, exprs_values = "logcounts", ncomponents = 200)
  sce_qc <- runTSNE(sce_qc, exprs_values = "logcounts", ntop = 500, ncomponents = 3)
  sce_qc <- runUMAP(sce_qc, exprs_values = "logcounts", ntop = 500, ncomponents = 3,
                    min_dist = 0.01, n_neighbors = 15, metric = "euclidean")
  
  # save the output 
  saveRDS(sce_qc, file = output_file_name)
}

normalize_sce(path_to_QCed_sce = args$path_to_QCed_sce, 
              seed = args$seed, 
              output_file_name = args$output_file_name, 
              HVG = args$HVG, 
              top_HVG = args$top_HVG)
