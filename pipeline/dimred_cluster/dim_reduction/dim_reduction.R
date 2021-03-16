perform_dim_reduction <- function(path_to_sce_norm, 
                                  HVG, 
                                  top_HVG, 
                                  top_PCs) {
  
  # load sce
  sce_norm <- readRDS(path_to_sce_norm)
  
  # compute the number of PCs to use 
  if (top_PCs == "computed") { # pipeline decides on the number of PCs to use 
    
    print("Computing PCA.")
    
    set.seed(1010)
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
    
    set.seed(101)
    sce_norm <- runPCA(sce_norm, exprs_values = "logcounts", ncomponents = PCs)
    
  } else { # compute PCA (and other dim reduction methods) without subsetting to top HVG
    
    set.seed(101)
    sce_norm <- runPCA(sce_norm, exprs_values = "logcounts", ncomponents = PCs)
    
  }
  
  # TSNE and UMAP 
  set.seed(101)
  print("Computing t-SNE.")
  sce_norm <- runTSNE(sce_norm, dimred = "PCA", exprs_values = "logcounts", ncomponents = 3, perplexity = 20)
  
  print("Computing UMAP.")
  sce_norm <- runUMAP(sce_norm, dimred = "PCA", exprs_values = "logcounts", ncomponents = 3, min_dist = 0.5, n_neighbors = 15, metric = "euclidean")
  
  print("Finished dimensionality reduction.")
  
  return(sce_norm)
  
} # end of function

