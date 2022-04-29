######### utils_dge.R
####
# This script serves as repo for DGE specific functions
# Author: Maxwell Douglas
# February 5th, 2020
#
# List of functions:
#
# fMarkersSampling()
# assembleMarkersByRank()
# assembleMarkersByP()
# plotDEGs()



# This function calculates differentially expressed genes between two clusters/samples of highly unequal sizes
# sce = a single cell experiment object
# pivot_cluster = The smaller sample or cluster
# target_group = The larger sample or cluster
# randomSubsets = The number of random samples to take of the larger cluster
# by doing repeated random sampling over the larger sample.
fMarkersSampling <- function(sce, pivot_cluster, target_group, randomSubsets) {
  
  ### Calculate the top DEGs for each random sample comparison
  idxs_list <- vector("list", randomSubsets)
  X <- 1:randomSubsets
  s_names <- sce$comparator
  idx_p <- which(s_names %in% c(pivot_cluster))
  idx_t <- which(s_names %in% c(target_group))
  
  U <-  list()
  for (i in X) { idxs_list[[i]] <- c(idx_p, sample(idx_t, size = length(idx_p))) }# grab random sample indices and pivot indices for loop
  
  cnt <- 1
  writeLines('\nRunning random sampling...')
  repeat {                                                                        # Do randomSubset iterations to find DEGs
    if(cnt > randomSubsets) {
      break
    }
    cat('.')
    idxs <- c(idx_p, sample(idx_t, size = length(idx_p)))
    x <- sce[rowData(sce)$qc_pass, idxs]                                          # Grab just the genes passing qc and cells of interest
    markers <- findMarkers(x, groups = x$comparator, log.p = TRUE)                        # Do pairwise differential expression (ANY)
    markers <- markers[[pivot_cluster]]
    markers_sig <- as.data.frame(markers) %>% rownames_to_column('gene_symbol') %>% dplyr::filter(log.FDR < -1.6)
    markers_sig <- dplyr::arrange(markers_sig, dplyr::desc(abs(markers_sig[,4]))) # order
    tt_top <- head(markers_sig, 2000)                                             # Identify top '2000' genes
    U[[cnt]] <- tt_top
    cnt <- cnt+1
  }
  
  ############### WHHHHHYYYYYYYY DOOOOESSSNNNTTTTTTT this WOOOOOOOOORRRRRKKKKKK????????
  # L <- lapply(ls, function(l){
  #   x <- sce[rowData(sce)$qc_pass, l[[1]]]
  #   marks <- findMarkers(x, groups = x$id, log.p = TRUE)
  # })
  ############### 
  
  return(U)                                                                       # Return a matrix of gene ranks in order of p-value by sample
}


### Assemble the top DEGs across all sampling comparisons
# Y = Matrix of gene ranks by FDR or p-value or LogFC for multiple comparisons
assembleMarkersByRank <- function(Y) {
  
  dt <- as.data.frame(table(unlist(Y)))
  lookup <- as.data.frame(integer(length(dt$Var1)))
  rownames(lookup) <- dt$Var1
  
  i <- 1
  while (i <= dim(Y)[2]) {                                                        # Count the ranks of each gene as it appears in the matrix
    n <- 1
    while (n <= dim(Y)[1]) {
      if (Y[n,i] %in% rownames(lookup)) {
        lookup[Y[n,i],] <- lookup[Y[n,i],] + n
      } else {
        break("No gene found")
      }
      n <- n+1
    }
    i <- i+1
  }
  
  test <- lookup
  lookup <- lookup %>% rownames_to_column('gene_symbol')
  lookup <- dplyr::arrange(lookup, lookup[,2])     # order
  # the_rank = data.frame(gene_rank = c(1:dim(dt)[1]))
  # test[,1] <- the_rank
  
  return(lookup)
}


### Assemble the top DEGs across all sampling comparisons using the median of adjusted p-values
# Y = List of dataframes containing the top DEGs for each sample comparison
assembleMarkersByP <- function(Y) {
  
  adjP <- lapply(Y, function(df) { df <- df %>% select(gene_symbol, log.FDR) } )
  fc <- lapply(Y, function(df) { df <- df %>% select(gene_symbol, names(Y[[1]])[5]) } )
  adjP <- adjP %>% purrr::reduce(full_join, by = c("gene_symbol" = "gene_symbol"))
  fc <- fc %>% purrr::reduce(full_join, by = c("gene_symbol" = "gene_symbol"))
  
  adjP_median <- apply(adjP[,2:ncol(adjP)], 1, median, na.rm = TRUE)
  fc_median <- apply(fc[,2:ncol(fc)], 1, median, na.rm = TRUE)
  
  degs <- data.frame(adjP$gene_symbol, adjP_median, fc_median)
  
  return(degs)
}

plotDEGs <- function(degs, clusters, ntop, reporter_gene) {
  
  names(degs) <- c('gene_symbol', 'Neg.Log.Adj.P', 'Log.FC')
  degs$Neg.Log.Adj.P <- -(degs[,2])
  degs <- degs[degs$gene_symbol != reporter_gene, ]                                          # there is no need to plot our reporter genes
  top <- head(degs, 50)
  p <- ggplot(degs, aes_string(x = names(degs)[3], y = names(degs)[2])) +
    geom_point(colour = 'gray') +
    geom_text_repel(data = top, aes(label = gene_symbol), colour = 'red', size = 5) +
    geom_point(data = top, aes_string(x = names(top)[3], y = names(top)[2]), colour = "blue") +
    labs(title = glue('{infected} vs. {control} top {ntop} DEGs from random sampling', infected = clusters[[1]], 
                      control = clusters[[2]], ntop = ntop),
         subtitle = paste('positive logFC means up-regulated in ', clusters[[1]]),
         caption = 'Adjusted P-vals and FCs are the median value from many random samples') + theme(plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}