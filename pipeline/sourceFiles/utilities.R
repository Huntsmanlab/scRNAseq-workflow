# SOURCE THIS FILE IN EVERY SINGLE SCRIPT YOU WRITE, OK? 

# LOAD PACKAGES ####
suppressPackageStartupMessages({
  
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(tidyverse)
  library(assertthat)
  library(BiocSingular)
  library(Matrix)
  library(argparse)
  library(data.table)
  library(cowplot)
  library(DropletUtils)
  library(ggplot2)
  library(gridExtra)
  library(styler)
  library(devtools)
  library(SC3)
  library(pheatmap)
  library(here)
  library(knitr)
  library(kableExtra)
  library(Seurat)
  library(batchelor)
  library(edgeR)
  library(limma)
  library(annotables)
  library(biomaRt)
  library(org.Hs.eg.db)
  library(EnsDb.Hsapiens.v75)
  library(AnnotationDbi)
  
})

#ADD THEME #### 
theme_amunzur <- theme(
  
  aspect.ratio = 1.0,
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line = element_line(size = 1),
  axis.line.x = element_line(color = "black", size = 1),
  axis.line.y = element_line(color = "black", size = 1),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(color = "black"),
  axis.title = element_text(color = "black"),
  axis.title.y = element_text(vjust = 0.2, size = 12),
  axis.title.x = element_text(vjust = 0.1, size = 12),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  legend.position = "none"
)


# before using the combine funtion, run this function first to turn our sces into seurat object. 
# run this function for each d

seuratClustering  <- function(sce, id) {
  
  # convert the sce into a seurat object 
  sobject <- as.Seurat(sce)
  
  # find variable features to reduce run time
  sobject <- FindVariableFeatures(sobject)
  
  # add to metadata so that we can keep track of the samples 
  sobject@meta.data[, "protocol"] <- paste(id)
  
  # we need to do this annoying thing where we do a bunch of renaming because seurat doesnt like capital letters 
  sobject@reductions$pca <- sobject@reductions$PCA
  sobject@reductions$tsne <- sobject@reductions$TSNE
  sobject@reductions$umap <- sobject@reductions$UMAP
  
  # compute the nearest neighbor graph
  sobject <-FindNeighbors(sobject)
  
  # lets find some clusters 
  sobject <- FindClusters(sobject, resolution = 0.15)
  
  return(sobject)
}


combine_sces2 <- function(seurat_list, ids_list) {
  
  # pick the first element of our seurat_list
  sobject1 <- seurat_list[[1]]
  
  # remove that object from the seurat_list
  seurat_list[[1]] <- NULL
  
  # convert the list to a vector, seurat needs this to be a vector 
  seurat_vector <- unlist(seurat_list)
  
  # convert the id list to a vector as well
  ids_vector <- unlist(ids_list)
  
  combined_seurat <-  merge(sobject1, y = seurat_vector, add.cell.ids = ids_vector, project = "protocol", data = TRUE)
  
  # convert this back to sce 
  combined_sce <- as.SingleCellExperiment(combined_seurat)
 
  return(combined_seurat) 
}

# # make seurat objects
# DH22control_seurat <- seuratClustering(DH22control, 'DH22control')
# DH22new_seurat <- seuratClustering(DH22new, 'DH22new')
# 
# DH21new_seurat <- seuratClustering(DH21new, 'DH21new')
# DH21control_seurat <- seuratClustering(DH21control, 'DH21control')
# 
# # combine them
# seurat_list <- list(DH22control_seurat, DH22new_seurat, DH21new_seurat, DH21control_seurat)
# ids_list <- list('DH22control', 'DH22new', 'DH21new', 'DH21control')
# 
# combined_sces <- combine_sces2(seurat_list, ids_list)
