library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "perform dim reduction in sce after normalization")

parser$add_argument('--path_to_sce_norm', metavar='DIRECTORY', type='character', help="Path normalized sce")

parser$add_argument('--output_file_name', metavar='FILE', type='character', help="Path to normalized sce")

parser$add_argument('--seed1', metavar='FILE', type='integer', help="seed for modelGeneVarByPoisson")

parser$add_argument('--seed2', metavar='FILE', type='integer', help="seed for UMAP, TNSE, PCA approximation")

parser$add_argument('--HVG', metavar = 'FILE', type = 'character', help = 'Do you want to subset to HVG?')

parser$add_argument('--top_HVG', metavar='FILE', type='integer', help="percentage for HVG")

parser$add_argument('--top_PCs', metavar='FILE', type='character', help="number PCs to used")


args <- parser$parse_args()

perform_dim_reduction <- function(path_to_sce_norm, 
                          output_file_name, 
                          seed1,
                          seed2,
                          HVG, 
                          top_HVG, 
                          top_PCs) {
  
  # load sce
  sce_norm <- readRDS(path_to_sce_norm)
  
  # compute the number of PCs to use 
  if (top_PCs == "computed") { # pipeline decides on the number of PCs to use 
    
    set.seed(seed1)
    dec.sce <- modelGeneVarByPoisson(sce_norm) # model gene variation, decompose it into technical and biological components
    
    # choose the HVG
    fraction_HGV <- top_HVG / 100
    top.sce <- getTopHVGs(dec.sce, prop = fraction_HGV)
    
    denoised.sce <- denoisePCA(sce_norm, technical = dec.sce)
    PCs <- ncol(reducedDim(denoised.sce))
    
  } else { # user decided how many PCs to use 
    
    top_PCs = as.numeric(top_PCs)
    PCs = top_PCs # number of PCs to use is given by the user in the top_PCs parameter
  }
  
  # deal with HVG
  if (HVG == "yes") { # user wants to subset the sce to HVG 
    
    # pick the HVG
    gene_var <- modelGeneVar(sce_norm)
    
    fraction_HGV <- top_HVG / 100
    chosen <- getTopHVGs(gene_var, prop = fraction_HGV)
    
    # subset to HVG
    sce_qc_hvg <- sce_norm[chosen,]
    
    # add the original count matrix with non HVG to our new sce
    altExp(sce_qc_hvg, "original") <- sce_norm
    sce_norm <- sce_qc_hvg
    
    set.seed(seed2)
    sce_norm <- runPCA(sce_norm, exprs_values = "logcounts", ncomponents = PCs)
    
  } else { # compute PCA (and other dim reduction methods) without subsetting to top HVG
    
    set.seed(seed2)
    sce_norm <- runPCA(sce_norm, exprs_values = "logcounts", ncomponents = PCs)
    
  }
  
  # TSNE and UMAP 
  set.seed(seed2)
  sce_norm <- runTSNE(sce_norm, dimred = "PCA", exprs_values = "logcounts", ncomponents = 3, perplexity = 20)
  sce_norm <- runUMAP(sce_norm, dimred = "PCA", exprs_values = "logcounts", ncomponents = 3, min_dist = 0.5, n_neighbors = 15, metric = "euclidean")
  
  # save the output 
  saveRDS(sce_norm, file = output_file_name)

} # end of function

perform_dim_reduction(path_to_sce_norm = args$path_to_sce_norm, 
              seed1 = args$seed1, 
              seed2 = args$seed2,
              output_file_name = args$output_file_name, 
              HVG = args$HVG, 
              top_HVG = args$top_HVG, 
              top_PCs = args$top_PCs)
