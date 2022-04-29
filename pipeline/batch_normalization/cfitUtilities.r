suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(cowplot)
  library(pheatmap)
  library(stringr)
  library(readxl)
  library(fgsea)
  library(tibble)
  library(ggplot2)
  library(scran)
  library(scater)
  library(glue)
  library(here)
  library(SingleR)
  library(cFIT)
  library(devtools)
  library(reshape2)
  library(destiny)
  library(Matrix)
  library(tidyverse)
  library(ggpointdensity)
  library(viridis)
  library(patchwork)
})

cFIT_plot_umap <- function(X = NULL, labels = NULL, pca = 50, n_components = 2, n_neighbors = 30, 
                      min_dist = 0.1, point.size = 0.3, alpha = 1, title = NULL, legend.name = "labels", 
                      cols = NULL, emb = NULL, seed = 0) {
  library(ggplot2)
  
  if (is.null(X) & is.null(emb)) {
    stop("data not provided!")
  }
  
  set.seed(seed)
  
  if (is.null(emb)) {
    if (!is.null(pca)) {
      if (pca > ncol(X)/2) {
        pca = NULL
      }
    }
    emb = uwot::umap(X, n_neighbors = n_neighbors, n_components = n_components, 
                     min_dist = min_dist, pca = pca)
  }
  
  if(n_components == 3) {
    df = data.frame(umap1 = emb[, 1], umap2 = emb[, 2], umap3 = emb[, 3], labels = if (!is.null(labels)) 
      labels else rep(0, nrow(X)))
    p = ggplot2::ggplot(df, aes(x = umap1, y = umap2, z = umap3)) + geom_point(col = "black", size = point.size, stroke = 0, shape = 16, alpha = alpha) + labs(x = "UMAP_1", y = "UMAP_2", title = title) + theme_light() + theme(plot.title = element_text(hjust = 0.5), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
  } else {
    df = data.frame(umap1 = emb[, 1], umap2 = emb[, 2], labels = if (!is.null(labels)) 
      labels else rep(0, nrow(X)))
    p = ggplot2::ggplot(df, aes(x = umap1, y = umap2)) + geom_point(col = "black", size = point.size, stroke = 0, shape = 16, alpha = alpha) + labs(x = "UMAP_1", y = "UMAP_2", title = title) + theme_light() + theme(plot.title = element_text(hjust = 0.5), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
  }
  
  if (!is.null(labels)) {
    if (is.null(legend.name)) {
      legend.name = "labels"
    }
    
    if (is.null(cols)) {
      cols = gg_color_hue(length(unique(labels)))
    }
    
    p = p + geom_point(aes(colour = labels), size = point.size, stroke = 0, shape = 16, 
                       alpha = alpha) + scale_color_manual(values = cols) + guides(col = guide_legend(ncol = 1, 
                                                                                                      title = legend.name, override.aes = list(size = 5)))
  }
  list(p = p, emb = emb)
}