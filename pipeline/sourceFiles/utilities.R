# SOURCE THIS FILE IN EVERY SINGLE SCRIPT 


# LOAD PACKAGES ####
# this is a list of all the libraries we need in our pipeline 
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
  library(org.Hs.eg.db)
  library(EnsDb.Hsapiens.v75)
  library(AnnotationDbi)
  library(gtools)
  library(formattable)
  library(DT)
  library(grid)
  library(ggplotify)
  library(EnhancedVolcano)
  library(ggrepel)
  
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

####################################################################################################
# list of functions
####################################################################################################
# seuratClustering()
# combine_sces2()
# make_upreg_table()
# make_downreg_table()
# make_pseudo_counts()
# find_reporter_genes()
# intersect_all()

# # from max's utilities folder: 
# combine_sces()
# gene_filter()
# plotDEG_scran()
# scran_batch_norm()
####################################################################################################

# this fucntion finds common genes in a given list of sces and subsets them to the common genes 
intersect_all <- function(sces){
  rownames_list <- lapply(sces, function(sce) rownames(sce)) # extract row names 
  universal <- Reduce(intersect, rownames_list) # find the common row names in the list 
  sces <- lapply(sces, function(sce) sce[universal, ]) # subset to common genes 
}


# before using the combine funtion, run this function first to turn our sces into seurat object. 
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

####################################################################################################

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

####################################################################################################

# this is to be used with DGE analysis with edge r. 
make_upreg_table <- function(qlf)  {
  
  up <- qlf$table[, c(1, 4)]
  
  # %>% 
  #   rownames_to_column() %>% # convert rownames to column to retain them in filter function 
  #   dplyr::filter(logFC > 0.5) %>% # only keep the genes with a logFC more than 0.5 
  #   dplyr::filter(PValue < 0.05) %>% # only keep the genes with p values less than 0.05 
  #   column_to_rownames() # convert the gene names back to rownames 
  
  
  # now we do the ranking. we will sort the logFC column from highest to lowest, gene at row1 (with the highest logFC) will be rank1, so the lowest rank. 
  # we will also sort the genes based on the p values, from lowest the highest. gene at row1 (with the lowest p value, most significant) will be rank1, so lowest rank. 
  # to decide how to report all the genes, we will use a combination of these two ranks. we will add the ranks and gene with the lowest added rank will be reported first. 
  
  logFC_sorted <- up[order(-up$logFC), ] # order the genes with large to small logFC
  logFC_sorted <- rownames_to_column(logFC_sorted)
  logFC_ranks <- c(1:nrow(logFC_sorted)) # assign the rank to genes
  logFC_df <- as.data.frame(cbind(logFC_sorted, logFC_ranks))
  logFC_df <- logFC_df[order(logFC_df$rowname),] # order alphabetically 
  
  PValue <- up[order(up$PValue), ] # order the genes with small to large p values 
  gene_list <- rownames(PValue) # extract genes
  PValue_ranks <- c(1:length(gene_list)) # assign the rank to genes
  PValue_df <- as.data.frame(cbind(gene_list, as.numeric(PValue_ranks))) # make a data frame 
  PValue_df <- PValue_df[order(PValue_df$gene_list),] # order alphabetically 
  
  combined <- cbind(logFC_df, PValue_df[, 2]) # combine both data frames 
  
  # turn both columns into numeric before doing the addition 
  combined[, 4] <- as.numeric(as.character(combined[, 4]))
  combined[, 5] <- as.numeric(as.character(combined[, 5]))
  
  # sum the ranks 
  combined$summed <- rowSums(combined[, 4:5]) 
  
  # reorder based on the summed ranks. remember that we are interested in genes with the lowest summed rank 
  combined <- combined[order(combined$summed), ]
  
  # make a table with this new order 
  interesting_genes <- combined %>% 
    dplyr::select(rowname, logFC, PValue, summed) %>% 
    dplyr::rename(summed_ranks = summed, gene_names = rowname) %>% 
    rownames_to_column() %>% 
    dplyr::select(gene_names, logFC, PValue, summed_ranks)
  
  return(interesting_genes)
  
}

####################################################################################################

make_downreg_table <- function(qlf)  {
  
  down <- qlf$table[, c(1, 4)]
  
  
  
  
  # %>%
  #   rownames_to_column() %>% # convert rownames to column to retain them in filter function
  #   dplyr::filter(logFC < -0.5) %>% # only keep the genes with a logFC more than 0.5
  #   dplyr::filter(PValue < 0.05) %>% # only keep the genes with p values less than 0.05
  #   column_to_rownames() # convert the gene names back to rownames
  
  
  # now we do the ranking. we will sort the logFC column from highest to lowest, gene at row1 (with the highest logFC) will be rank1, so the lowest rank. 
  # we will also sort the genes based on the p values, from lowest the highest. gene at row1 (with the lowest p value, most significant) will be rank1, so lowest rank. 
  # to decide how to report all the genes, we will use a combination of these two ranks. we will add the ranks and gene with the lowest added rank will be reported first. 
  
  logFC_sorted <- down[order(down$logFC), ] # order the genes with large to small logFC
  logFC_sorted <- rownames_to_column(logFC_sorted)
  logFC_ranks <- c(1:nrow(logFC_sorted)) # assign the rank to genes
  logFC_df <- as.data.frame(cbind(logFC_sorted, logFC_ranks))
  logFC_df <- logFC_df[order(logFC_df$rowname),] # order alphabetically 
  
  PValue <- down[order(down$PValue), ] # order the genes with small to large p values 
  gene_list <- rownames(PValue) # extract genes
  PValue_ranks <- c(1:length(gene_list)) # assign the rank to genes
  PValue_df <- as.data.frame(cbind(gene_list, as.numeric(PValue_ranks))) # make a data frame 
  PValue_df <- PValue_df[order(PValue_df$gene_list),] # order alphabetically 
  
  combined <- cbind(logFC_df, PValue_df[, 2]) # combine both data frames 
  
  # turn both columns into numeric before doing the addition 
  combined[, 4] <- as.numeric(as.character(combined[, 4]))
  combined[, 5] <- as.numeric(as.character(combined[, 5]))
  
  # sum the ranks 
  combined$summed <- rowSums(combined[, 4:5]) 
  
  # reorder based on the summed ranks. remember that we are interested in genes with the lowest summed rank 
  combined <- combined[order(combined$summed), ]
  
  # make a table with this new order 
  interesting_genes <- combined %>% 
    dplyr::select(rowname, logFC, PValue, summed) %>% 
    dplyr::rename(summed_ranks = summed, gene_names = rowname) %>% 
    rownames_to_column() %>% 
    dplyr::select(gene_names, logFC, PValue, summed_ranks)
  
  return(interesting_genes)
  
}

####################################################################################################

# make pseudo counts 
make_pseudo_counts <- function(sce_list) {
  
  i <-  1 
  pseudoCountsList <-  list() # make an empty list to append the results of the while loop 
  
  # make the pseudo bulk counts: 
  while (i <= length(sces)){
    
    sce <- sces[[i]] # pick the sce from the sce list 
    combo <- combos[[i]] # pick the related element for the index we are at. combos is a guide, tells us how to do the summing. 
    sum <- sumCountsAcrossCells(counts(sce), combo) # sum across cluster. this means add all the expression values for each gene for all cells in the cluster 
    
    colnames(sum) <- as.numeric(colnames(sum)) # convert the column names to integers
    
    sum <- as.data.frame(sum) # now convert it from double to data frame. 
    
    colnames(sum) <- paste(colnames(sum), paste(dh_organoid_ids[[i]]), sep = '_') # add the id of the sample next to the cluster number
    
    pseudoCountsList <- append(pseudoCountsList, list(sum)) # add your results to the list and update it at each iteration. this list has the new counts that we are interested in. 
    
    i <- i + 1 # go to the next index 
    
  } # end of the while loop
  
  return(pseudoCountsList)
  
} # end of function 

####################################################################################################


find_reporter_genes <- function(reporter_name, sce) {
  
  index <- grep(reporter_name, rownames(counts(sce)), ignore.case=TRUE) # find the reporter gene index  
  reporter_counts <- as.data.frame((counts(sce)[index, ])) # extract counts for the reporter gene 
  
  reporter_num <- apply(reporter_counts, 2, function(c) sum(c != 0)) # find the number of cells that express the reporter gene 
  
  return(reporter_num)
  
}


####################################################################################################


# This function creates and returns a plot of differentially expressed genes between two samples
plotDEG_scran <- function(ids, markers, ntop) {
  
  # Massage output from DE calculations into useful format
  markers <- markers[[ids[1]]]                                                        # A +ve logFC = up-regulated in tranduced cells
  markers_all <- as.data.frame(markers) %>% rownames_to_column('gene_symbol')
  
  if ("FDR" %in% colnames(markers_all)) {                                             # Case when p = p.value
    markers_sig <- as.data.frame(markers) %>% rownames_to_column('gene_symbol') %>% dplyr::filter(FDR <= 0.05)
    tooLow <- function(x) (ifelse(x < 1e-315, 1e-315, x))
    markers_sig <- markers_sig %>% mutate_at(c("p.value", "FDR"), tooLow)
    # Order the genes for plotting
    markers_sig <- dplyr::arrange(markers_sig, dplyr::desc(abs(markers_sig[,4])))     # order
    t_mark <- markers_sig                                                             # make copy
    t_mark$FDR <- -log10(t_mark[,3])                                                  # log transform
    tt_top <- head(t_mark, ntop)                                                      # Identify top 'n' genes
    
  } else {                                                                            # Case when p = log.p.value
    markers_sig <- as.data.frame(markers) %>% rownames_to_column('gene_symbol') %>% dplyr::filter(log.FDR < -1.6)
    # Order the genes for plotting
    markers_sig <- dplyr::arrange(markers_sig, dplyr::desc(abs(markers_sig[,4])))     # order
    t_mark <- markers_sig                                                             # make copy
    t_mark$log.FDR <- -(t_mark[,3])
    tt_top <- head(t_mark, ntop)                                                      # Identify top 'n' genes
  }
  
  
  # Plot the differentially expressed genes
  p1 <- ggplot(t_mark, aes_string(x = names(t_mark)[4], y = names(t_mark)[3])) +
    geom_point(colour = 'gray') +
    geom_text_repel(data = tt_top, aes(label = gene_symbol), colour = 'red', size = 5) +
    geom_point(data = tt_top, aes_string(x = names(t_mark)[4], y = names(t_mark)[3]), colour = "blue") +
    labs(title = glue('{infected}  vs. {control} top {ntop} DEGs', infected = ids[[1]], 
                      control = substr(names(t_mark)[4], 7, nchar(names(t_mark)[4])), ntop = ntop),
         subtitle = 'positive logFC means up-regulated in transduced sample') + theme(plot.subtitle = element_text(hjust = 0.5))
  
  markers_sig <- dplyr::arrange(markers_sig, dplyr::desc(abs(markers_sig[,5])))       # order genes for return matrix
  t_mark <- dplyr::arrange(t_mark, dplyr::desc(abs(t_mark[,5])))
  tt_top <- head(t_mark, ntop)
  
  p2 <- ggplot(t_mark, aes_string(x = names(t_mark)[5], y = names(t_mark)[3])) +
    geom_point(colour = 'gray') +
    geom_text_repel(data = tt_top, aes(label = gene_symbol), colour = 'red', size = 5) +
    geom_point(data = tt_top, aes_string(x = names(t_mark)[5], y = names(t_mark)[3]), colour = "blue") +
    labs(title = glue('{infected}  vs. {control} top {ntop} DEGs', infected = ids[[1]], 
                      control = substr(names(t_mark)[5], 7, nchar(names(t_mark)[5])), ntop = ntop),
         subtitle = 'positive logFC means up-regulated in first sample') + theme(plot.subtitle = element_text(hjust = 0.5))
  
  p <- plot_grid(p1, p2, labels = "AUTO")
  
  outputs <- list(markers_all = markers_all, markers_sig = markers_sig, tt_top = tt_top, plot = p)
}

#####
# Do batch-correction on several batches or samples
scran_batch_norm <- function(sce_list) {
  
  # Make sure each sce has the same number of rows (we will do so by subsetting to the common rows)
  row_list <- lapply(sce_list, function(sce) rownames(sce))
  common_rows <- Reduce(intersect, row_list)
  if(length(common_rows) == 0){
    stop('SingleCellExperiment objects don\'t share any common genes, cannot combine cells together.')
  }
  # extract common genes
  sce_list <- lapply(sce_list, function(sce, keep_rows){
    sce <- sce[keep_rows, ]
  }, keep_rows = common_rows )                                                  
  
  sce_list <- do.call(batchelor::multiBatchNorm, c(sce_list, assay.type="counts", min.mean = 0.1))
  
  sce_list <- lapply(sce_list, function(sce){
    for(n in reducedDimNames(sce)){
      reducedDim(sce, n) <- NULL
    }
    return(sce)
  })
  return(sce_list)
}

####################################################################################################
# Filter genes out based on their detection rate across all pseudocells (or 'libraries' in the common scRNA parlance)
# Good default is rm genes whose total count < 100 and detected in < 5% of cells, remove mito, ribo genes
# Example: gene_filter(sce, gene_min_counts = 100, gene_min_detection_rate = 5, reporters = NULL)
gene_filter <- function(sce, gene_min_counts, gene_min_detection_rate, reporters = NULL){
  
  count_drop <- rowSums(counts(sce)) < gene_min_counts
  rate_drop <- (rowSums(counts(sce) > 0) / dim(sce)[2] * 100) < gene_min_detection_rate
  mt_genes <- grepl("^MT-", rowData(sce)$Symbol)
  ribo_genes <- grepl("^RP[LS]", rowData(sce)$Symbol)
  
  gene_drop <- (count_drop | rate_drop | mt_genes | ribo_genes) 
  
  if(!is.null(reporters)) {
    report_genes <- rownames(sce) %in% reporters
    gene_drop <- (gene_drop | report_genes )
  }
  
  rowData(sce)$qc_pass <- !gene_drop
  
  return(sce)
  
}

####################################################################################################
# Combine multiple SCEs together
combine_sces <- function(..., sce_list = NULL, prefix_col_name = NULL, suffix_colnames = NULL, prefix_cell_name = FALSE, gene_meta_colnames = NULL){
  
  if(is.null(sce_list)){
    sce_list <- list(...)
  }
  
  if(length(sce_list) == 1){
    sces <- sce_list[[1]]
    
  } else{ # We have multiple sces! Let's do it
    
    cell_meta_colnames <- lapply(sce_list, function(sce) colnames(colData(sce)))
    common_cell_meta_cols <- Reduce(intersect, cell_meta_colnames)
    if(length(cell_meta_colnames) == 0){
      stop('SingleCellExperiment objects don\'t share any common cell_meta, cannot combine cells together.')
    }
    
    # cbind for the SingleCellExperiment package doesn't like conflicting meta data for rows. 
    # Logicals are fine. (I.E. qc_pass True/False)
    row_meta <- lapply(sce_list, function(sce) colnames(elementMetadata(rowRanges(sce))))
    gene_meta_colnames <- Reduce(intersect, row_meta)
    
    
    sce_list <- lapply(sce_list, function(sce, common_cell_meta_cols, prefix_col_name, suffix_colnames, prefix_cell_name, gene_meta_colnames){
      
      colData(sce) <- colData(sce)[, common_cell_meta_cols]                     # extract common cell meta
      metadata(sce) <- metadata(sce_list[[1]])                                  # make general sample metadata consistent
      
      # add id to values of selected cell metadata, e.g., 0 -> DH4-0 for cluster name
      if( (!is.null(prefix_col_name)) | (!is.null(suffix_colnames)) ){
        for(colname in suffix_colnames){
          colData(sce)[, colname] <- paste(colData(sce)[, prefix_col_name], colData(sce)[, colname], sep = '.')
        }
      }
      
      # add id to cell names because different datasets can contain the same cell name
      if(prefix_cell_name & !is.null(prefix_col_name)){
        colnames(sce) <- paste(colData(sce)[, prefix_col_name], colnames(sce), sep = '.')
      }
      
      # cbind needs a consistent number of dimensions for the projections, we will truncate to the minimum (2 for each)
      if (!is.null(sce_list[[1]]@int_colData@listData[["reducedDims"]])) {
        sce@int_colData@listData[["reducedDims"]]@listData <- lapply(sce@int_colData@listData[["reducedDims"]]@listData, function(list) {list <- list[,1:2]})
      }
      # Gene metadata must be the same
      # To do that we are going to use the values from one, for all. Thus overwriting the conflicting metadata...
      if(!is.null(gene_meta_colnames)){
        rowData(sce)[,gene_meta_colnames] <- rowData(sce_list[[1]])[, gene_meta_colnames]
      }
      
      ## Alternatively we could remove the conflicting gene metadata using the following
      # Note that this will also remove logicals associated with each gene
      # if(!is.null(gene_meta_colnames)){
      #   to_keeps <- intersect(gene_meta_colnames, colnames(rowData(sce)) )
      #   if(length(to_keeps) == 0){
      #     stop('No columns in rowData(x) has the same names as column names provided in gene_meta_colnames')
      #   }
      #   rowData(sce) <- rowData(sce)[, to_keeps]
      # }
      
      return(sce)
      
    }, common_cell_meta_cols = common_cell_meta_cols, prefix_col_name = prefix_col_name, suffix_colnames = suffix_colnames,
    prefix_cell_name = prefix_cell_name, gene_meta_colnames = gene_meta_colnames)
    
    # ISSUES WITH cbind:
    # 1) The docs say that Row and Experimental metadata will only be taken from the first element in the list 
    # but the function fails when sces contain any difference in gene metadata. (Other than logical masks)
    # 2) Fails if reduced dimension projections are of inconsistent dimensionality (you can only keep 2 components 
    # of the pca in all samples if just one had 2)
    # 3) Fails if colData contains different column names (this is ok)
    sces <- sce_list
    while (length(sce_list) >= 2) {
      sce_list[1] <- do.call(SingleCellExperiment::cbind, sce_list[1:2])
      sce_list[2] <- NULL 
    }
    
    
    # reducedDim(sces, "no_correction_PCA") <- reducedDim(sces, "PCA")
    # reducedDim(sces, "no_correction_UMAP") <- reducedDim(sces, "UMAP")
    # reducedDim(sces, "no_correction_TSNE") <- reducedDim(sces, "TSNE")
    # reducedDim(sces, "PCA") <- NULL
    # reducedDim(sces, "UMAP") <- NULL
    # reducedDim(sces, "TSNE") <- NULL
    
  }
  
  return(sce_list[[1]])
  
}