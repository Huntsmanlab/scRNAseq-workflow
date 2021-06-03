library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scran)
})


## Assign cell cycle phases to clustered sce and save the new sce object as cell_cycle.rds.
## The function we use is cyclone from scran.

assign_cell_cycle <- function(path_to_sce, species, output_path, corrected_path){
  
  # read clustered sce
  sce <- readRDS(path_to_sce)
  
  # read marker genes for cell cycle
  if (species=="human") {
    marker.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
  } else { marker.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))}
  
  # assign cell cycle phase to each cell; save the assignment information in the cell_cycle column
  print("assigning cell cycle phases.. this could take a while depending on cell numbers")
  assignments <- cyclone(sce, marker.pairs, gene.names=rowData(sce)$ID)
  sce$cell_cycle <- assignments[[1]]
  
  
  # for each gene, get the percentage of variance explained by cell cycle; keep those genes with %variance explained less than 5%
  diff <- getVarianceExplained(sce, DataFrame(sce$cell_cycle))
  discard <- diff > 5

  # label identified cell cycle genes as TRUE in the cell_cycle_gene column of sce rowData:
  # we can use this column to filter out cell cycle related genes when performing DGE analysis
  rowData(sce)$cell_cycle_gene <- as.logical(discard)
  rowData(sce)$cell_cycle_gene[is.na(rowData(sce)$cell_cycle_gene)] <- FALSE
  
  # save sce with cell cycle gene label information
  saveRDS(sce, file=output_path)
  
  
  # correct the cell cycle effect: treat cell cycle as a batch effect
  corrected <- multiBatchNorm(sce, batch=sce$cell_cycle)
  
  # save sce with batch (cell cycle) corrected log-counts
  saveRDS(corrected, file=corrected_path)

  out <- list(cycle_gene_label = sce, cycle_corrected = corrected)
  
  return(out)
  
}
  





