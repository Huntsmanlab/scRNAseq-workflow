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
  #library(annotables)
  library(org.Hs.eg.db)
  library(EnsDb.Hsapiens.v75)
  library(AnnotationDbi)
  library(gtools)
  library(formattable)
  library(DT)
  library(grid)
  library(ggplotify)
  #library(EnhancedVolcano)
  library(ggrepel)
  library(glue)
  library(fgsea)
  library(splatter)
  library(viridis)
  library(ggthemes)
  library(destiny)
  library(ggbeeswarm)
  library(plotly)
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
# make_pseudo_counts_seurat()
# find_and_bind() ----- find common rows in two or more data frames, subset and cbind the common rows 
# find_reporter_genes()
# intersect_all()
# seurat_integrate()

# # from max's utilities folder: 
# combine_sces()
# gene_filter()
# plotDEG_scran()
# scran_batch_norm()
####################################################################################################
# make_tsne_plot <- function(sce, 
#                            color, 
#                            dimensions, 
#                            title) {
#   
#   # extract coordinates and save them in a list 
#   i <- 1
#   dimensions_list <- list()
#   
#   while (i <= dimensions){
#     dim <- list(as.numeric(sce@int_colData@listData[["reducedDims"]]$TSNE[, i]))
#     dimensions_list <- append(dimensions_list, dim) 
#     
#     i <- i + 1 
#     
#   }
#   
#   # make a df from the coordinates 
#   df <- do.call(data.frame, dimensions_list)
#   
#   # do some renaming 
#   if (dimensions == 3){
#     names(df) <- c("t-SNE1", "t-SNE2", "t-SNE3")
#   } else{
#     names(df) <- c("t-SNE1", "t-SNE2")
#   }
#   
#   # now add another column to df based on what you want to color 
#   if (color == 'cluster') {
#     df$cluster <- sce$cluster
#   } else if (color == 'cell_type'){
#     df$cell_type <- sce$cell_type
#   }
#   
#   # make the plot here
#   if (dimensions == 3) {
#     
#     dimred_plot <- plot_ly(data = df,
#                            x = ~df[, 1], y = ~df[, 2], z = ~df[, 1]
#                            opacity = 1,
#                            color = ~cell_type,
#                            type = "scatter",
#                            mode = "markers",
#                            marker = list(size = 5)) %>% 
#       
#       add_markers()
#     
#     cell_type_plot <- cell_type_plot %>% 
#       layout(title = 'cell types in the sample shown in t-SNE plot')
#     
#     
#   }
#   
#   
#   return(dimred_plot)
#   
# }
# 
# 
# 
# 



# CELL ASSIGN COLORS
library(RColorBrewer)
n <- 20
qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
most_distinct_color_palette <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# 14 colors, one color for each cell type. cell types given in alphabetical order 
master_cell_types_old <- c("B cells",
                       "Cytotoxic T cells",
                       "Endometrial stem cells",
                       "Endothelial cells",            
                       "Epithelial cells",
                       "Epithelial ciliated cells",    
                       "Mesenchymal cells",            
                       "Mesenchymal stem cells",
                       "Monocyte/Macrophage",
                       "Myofibroblast",
                       "other",                        
                       "Plasma cells",
                       "Proliferative cells",
                       "T cells",
                       "Vascular smooth muscle cells", 
                       "Cluster2_like")

master_cell_types <- c("endothelial_cells",
                       "epithelial_cells",
                       "lymphocyte_cells",
                       "macrophage_cells",
                       "other",
                       "plasma_cells",
                       "stromal_fibroblast_cells")


master_cell_types_GC <- c("endothelial_cells",
                          "epithelial_cells",
                          "granulosa_cells",
                          "lymphocyte_cells",
                          "macrophage_cells",
                          "mesothelial_cells",
                          "other",
                          "plasma_cells",
                          "stromal_fibroblast_cells",
                          "theca_cells")

master_color_palette <- c("gold1", 
                          "lightskyblue",
                          "darkslateblue",
                          'thistle4',
                          'mediumorchid2',
                          "indianred1",
                          "plum2",
                          "springgreen",
                          'royalblue1',
                          'chartreuse1',
                          "darkmagenta")

# master_color_palette <- colorRampPalette(brewer.pal(8, "Set1"))(20)

# to visualize these colors 
# pie(rep(1, 15), col = master_color_palette)

# this function plots cell assign results in a dim reduction plot. 
visualize_cellassign <- function(seurat_object, reduction_type, group_by, master_color_palette, master_cell_types){
  
  # find cell types present in the seurat object, sort them alphabetically
  unique_cell_types <- as.vector(sort(unique(seurat_object$cell_type)))
  
  # find common cell types between our lists and the master cell type list 
  common <- intersect(unique_cell_types, master_cell_types)
  
  # find the indices where we have common these cell types in the master cell type list 
  idx <- match(common, master_cell_types)
  
  # pick colors that correspond to the indices
  colors_chosen <- master_color_palette[idx]
  
  # make cell assign plot
  plot <- DimPlot(object = seurat_object,
                  dims = c(1, 2),
                  reduction = reduction_type,
                  group.by = group_by,
                  pt.size = 1.5) +
    scale_color_manual(values = colors_chosen)
  
  return(plot)
  
} # end of function

# this fucntion finds common genes in a given list of sces and subsets them to the common genes 
intersect_all <- function(sces){
  rownames_list <- lapply(sces, function(sce) rownames(sce)) # extract row names 
  universal <- Reduce(intersect, rownames_list) # find the common row names in the list 
  sces <- lapply(sces, function(sce) sce[universal, ]) # subset to common genes
  
  return(sces)
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

make_pseudo_counts_seurat <- function(sobject) {
  
  new_names <- strsplit(colnames(sobject), '.', fixed = TRUE)
  new_names <- lapply(new_names, function(element) as.list(element)[[1]])
  new_names <- c(unlist(new_names))
  
  sobject_df <- as.data.frame(GetAssayData(sobject))
  colnames(sobject_df) <- new_names
  
  # here we sum the gene counts for the rows that have the same column name. this allows us to sum all the counts for a specific sample 
  sobject_df <- sapply(unique(colnames(sobject_df)), function(x) rowSums(sobject_df[, colnames(sobject_df) == x, drop=FALSE]))
  
  sobject_df <- as.data.frame(sobject_df)
  
  return(sobject_df)
}

####################################################################################################

make_pseudo_counts_sce <- function(sce_list) {
  
  # subset all sces to common genes 
  
  df_list <- do.call(counts, sce_list)
  df_sum <- rowSums(df)
  
  return(df_sum)
  
}


####################################################################################################

# assuming that we already subsetted sces to common genes before, so make sure to run intersect_all
find_and_bind_multiple <- function(sce_list) {
  
  # extract counts - make a df 
  df_list <- lapply(sce_list, function(sce) as.data.frame(as.matrix(counts(sce))))
  
  # now sum the rows 
  df_list <- lapply(df_list, function(df) as.data.frame(rowSums(df)))
  
  # turn rownames to columns for all data frames
  df_list <- lapply(df_list, function(df) df %>% rownames_to_column())
  
  # do an inner join by rowname
  combined <- Reduce(function(x, y) inner_join(x, y, by = "rowname"), df_list)
  
  # col back to rowname 
  combined <- combined %>% column_to_rownames()
  
  # rename the columns 
  combined_sce <- do.call(combine_sces, sce_list)
  names <- unique(combined_sce$id)
  colnames(combined) <- names
  
  return(combined)
  
}

####################################################################################################

find_and_bind <- function(df_list){
  
  # find common genes
  universal <- intersect(rownames(df_list[[1]]), rownames(df_list[[2]]))
  
  # subset to common genes 
  df_list[[1]] <- df_list[[1]][universal, ]
  df_list[[2]] <- df_list[[2]][universal, ]
  
  # order genes alphabetically so that samples will have the same order of genes 
  df_list[[1]] <- df_list[[1]][ order(row.names(df_list[[1]])), ]
  df_list[[2]] <- df_list[[2]][ order(row.names(df_list[[2]])), ]
  
  # now cbind: 
  combined <- cbind(df_list[[1]], df_list[[2]])
  
  # remove negative counts
  combined <- combined[apply(combined, 1, function(combined) all(combined >= 0)), ]
  
  return(combined)
  
}

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

fMarkersSampling <- function(sce, pivot_cluster, target_group, randomSubsets) {
  
  ### Calculate the top DEGs for each random sample comparison
  idxs_list <- vector("list", randomSubsets)
  X <- 1:randomSubsets
  s_names <- sce$id
  idx_p <- which(s_names %in% c(pivot_cluster))
  idx_t <- which(s_names %in% c(target_group))
  
  U <-  list()
  for (i in X) { idxs_list[[i]] <- c(idx_p, sample(idx_t, size = length(idx_p))) }# grab random sample indices and pivot indices for loop
  
  cnt <- 1
  repeat {                                                                        # Do randomSubset iterations to find DEGs
    if(cnt > randomSubsets) {
      break
    }
    
    idxs <- c(idx_p, sample(idx_t, size = length(idx_p)))
    x <- sce[rowData(sce)[[3]], idxs]                                             # Grab just the genes passing qc and cells of interest
    markers <- findMarkers(x, groups = x$id, log.p = TRUE)                        # Do pairwise differential expression (ANY)
    markers <- markers[[pivot_cluster]]
    markers_sig <- as.data.frame(markers) %>% rownames_to_column('gene_symbol') %>% dplyr::filter(log.FDR < -1.6)
    markers_sig <- dplyr::arrange(markers_sig, dplyr::desc(abs(markers_sig[,4]))) # order
    tt_top <- head(markers_sig, 2000)                                             # Identify top '2000' genes
    
    U[[cnt]] <- tt_top
    cnt <- cnt+1
  }
  
  return(U)                                                                       # Return a matrix of gene ranks in order of p-value by sample
}

########################################################################################################################

# integrate samples with seurat integration; returns integrated sce and combined sce without batch normalization

seurat_integrate <- function(path_to_sce_clus, 
                             ids_integration) {
  
  id.list <- ids_integration
  
  id.list <- strsplit(id.list, "-")[[1]] # split by "-" SO ID must NOT CONTAIN "-"
  id.list <- as.list(c(unlist(list(id.list))))
  id.list <- do.call(list, id.list)  
  
  id.orig <- id.list
  
  # lets start! 
  # load the data
  sces <- lapply(path_to_sce_clus, function(path) readRDS(path))
  
  # subset to common genes across a group of sces 
  sces <- intersect_all(sces)
  
  # add the id number in front of cell barcodes
  repeated <- lapply(sces, function(sce) dim(sce)[[2]]) # extract cell numbers from each sce
  id.list <- mapply(rep, id.list, repeated) # repeat the cell ids as many times as the cell number 
  
  # sometimes different samples use same barcodes, we will append the id name in front of the barcodes to avoid that 
  # we use a period to separate the id and the barcode because some ids use '_' already 
  for (i in 1:length(sces)) {
    colnames(sces[[i]]) <- paste(id.list[[i]], colnames(sces[[i]]), sep = '.')
  }
  
  # convert them all to seurat objects 
  # note that we aren't using the normalized data since we will use seurat's normalization method 
  # counts(sce) gives raw counts
  seurats <- lapply(sces, function(sce) CreateSeuratObject(counts = counts(sce), min.cells = 3, min.features = 200)) 
  
  # normalize each sample using seurat's methods
  seurats <- lapply(seurats, function(seurat) NormalizeData(seurat, verbose = FALSE))
  
  # find variable features for each of our seurat objects 
  seurats <- lapply(seurats, function(seurat) FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000))
  
  # select the number of features to integrate 
  id.features <- SelectIntegrationFeatures(object.list = seurats, nfeatures = nrow(seurats[[1]]))
  
  # make a reference list for the next step  
  anchors <- FindIntegrationAnchors(object.list = seurats, dims = 1:30, k.filter = NA)
  
  # find the common genes 
  total.genes <- lapply(seurats, function(seurat) rownames(seurat@assays$RNA@counts))
  common.genes <- Reduce(f = intersect, x = total.genes)
  
  # now merge the data, passing the common genes help us only integrate those and disregard the other genes 
  integrated <- IntegrateData(anchorset = anchors, dims = 1:30, features.to.integrate = common.genes)
  
  # now we will add the sample ids each sample is associated with 
  sample_names <- strsplit(colnames(integrated), '.', fixed = TRUE)
  sample_names <- lapply(sample_names, function(element) as.list(element))
  sample_names <- lapply(sample_names, function(element) element[[1]])
  
  sample_names <- do.call(rbind, sample_names)
  rownames(sample_names) <- colnames(integrated)
  
  # add this to the metadata 
  integrated <- AddMetaData(integrated, sample_names, col.name = 'id')
  
  # the default assay is the new one 
  DefaultAssay(integrated) <- "integrated"
  
  # then some standard workflow for visualization and dim reduction 
  integrated <- ScaleData(integrated, verbose = FALSE)
  set.seed(1998)
  integrated <- RunPCA(integrated, verbose = FALSE)
  integrated <- RunTSNE(integrated, dims = 1:30, dim.embed = 3, seed.use = 300, perplexity = 20)
  integrated <- RunUMAP(integrated, dims = 1:30, n.components = 3L, seed.use = 1000)
  
  # clustering
  integrated <- FindNeighbors(integrated, dims = 1:10)
  integrated <- FindClusters(integrated, resolution = 0.5)
  
  
  # UNCORRECTED DATA
  # here we do the combined but uncorrected sample, with the batch effects
  # load the normalized data 
  sces_norm <- lapply(path_to_sce_clus, function(path) readRDS(path))
  
  sces_norm <- intersect_all(sces_norm) # subset to common genes 
  combined <- do.call(combine_sces, sces_norm) # combine the objects by doing a simple r bind 
  
  # some dim reduction 
  set.seed(1564)
  combined <- runPCA(combined)
  combined <- runTSNE(combined)
  combined <- runUMAP(combined)
  
  out <- list(corrected = integrated, uncorrected = combined)
  return(out)
  
} # end of function 




