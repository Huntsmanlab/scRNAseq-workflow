library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scran)
})


## Assign cell cycle phases to clustered sce and save the new sce object as cell_cycle.rds.
## The function we use is cyclone from scran.

assign_cell_cycle <- function(path_to_sce, cell_cycle_csv){
  
  # read clustered sce
  sce <- readRDS(path_to_sce)
  
  # read marker genes for cell cycle
  if(grepl("^ENSMUS", sce@rowRanges@elementMetadata@listData[["ID"]][1])){
    # mouse sample
    marker.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  } else {
    # human sample
    marker.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
  }
  
  # assign cell cycle phase to each cell; save the assignment information in the cell_cycle column
  print("assigning cell cycle phases.. this could take a while depending on cell numbers")
  assignments <- cyclone(sce, marker.pairs, gene.names=rowData(sce)$ID)
  sce$cell_cycle <- assignments[[1]]
  
  # put barcodes and cell cycle phases into dataframe
  cell_cycle <- as.vector(assignments[[1]])
  barcodes <- as.vector(colnames(sce)) 
  df <- as.data.frame(cbind(cell_cycle, barcodes))
  names(df) <- c("cell_cycle", "barcodes")
  
  # save cell cycle information into a csv file
  dir.create(dirname(cell_cycle_csv))
  write_csv(df, cell_cycle_csv) 
  
  
  # for each gene, get the percentage of variance explained by cell cycle; keep those genes with %variance explained less than 5%
  diff <- getVarianceExplained(sce, DataFrame(sce$cell_cycle))
  discard <- diff > 5

  # label identified cell cycle genes as TRUE in the cell_cycle_gene column of sce rowData:
  # we can use this column to filter out cell cycle related genes when performing DGE analysis
  rowData(sce)$cell_cycle_gene <- as.logical(discard)
  rowData(sce)$cell_cycle_gene[is.na(rowData(sce)$cell_cycle_gene)] <- FALSE
  
  
  # correct the cell cycle effect: treat cell cycle as a batch effect
  # corrected <- multiBatchNorm(sce, batch=sce$cell_cycle)
  
  return(sce)
  
}
  





