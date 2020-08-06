# This script is for integrating multiple data sets. integrating allows us to combine data sets and remove batch effects at the same time 
# we may need to integrate data in various situations, so this script can be used in multiple types of analysis, such as: 
# running a clustering analysis on the combined sample 
# run dim reduction plots on combined plots 
# run DGE analysis across samples, we will need to remove batch effects for that. we can integrate data sets first here, then split if needed. 
# data here is saved into the data/integrated folder/<name of integrated samples>
# this is the version integrated into snakemake 

# load a few necessary things 
library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "integrate multiple datasets to remove batch effects")

parser$add_argument('--path_to_sce_cas', nargs = '+', metavar = 'DIRECTORY', type = 'character',
                    help = "Path to sce cell assign analysis")

parser$add_argument('--path_to_combined_sce', metavar = 'DIRECTORY', type = 'character',
                    help= "Path to sce cell assign analysis")

parser$add_argument('--integrated_object_name', nargs = 1, metavar = 'FILE', type = 'character',
                    help = "Path to integrated samples, doesnt have batch effects")

parser$add_argument('--ids_integration', metavar = 'FILE', type = 'character',
                    help = "the name of the samples you want to integrate")

args <- parser$parse_args()


############################################################################################################################


whatagreatfunction <- function(path_to_sce_cas, 
                               path_to_combined_sce, 
                               integrated_object_name, 
                               ids_integration) {
  
  id.list <-ids_integration
  
  id.list <- strsplit(id.list, "-")[[1]] # split by "-" SO ID must NOT CONTAIN "-"
  id.list <- as.list(c(unlist(list(id.list))))
  id.list <- do.call(list, id.list)  
  
  id.orig <- id.list
  
  # lets start! 
  # load the data, note that we aren't using the normalized data since we will use seurat's normalization method 
  sces <- lapply(path_to_sce_cas, function(path) readRDS(path))
  
  # subset to common genes across a group of sces 
  intersect_all <- function(sces){
    rownames_list <- lapply(sces, function(sce) rownames(sce)) # extract row names 
    universal <- Reduce(intersect, rownames_list) # find the common row names in the list 
    sces <- lapply(sces, function(sce) sce[universal, ]) # subset to common genes 
  }
  
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
  seurats <- lapply(sces, function(sce) CreateSeuratObject(counts = counts(sce), min.cells = 3, min.features = 200)) 
  
  # normalize each sample using seurat's methods
  seurats <- lapply(seurats, function(seurat) NormalizeData(seurat, verbose = FALSE))
  
  # find variable features for each of our seurat objects 
  seurats <- lapply(seurats, function(seurat) FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000))
  
  # select the number of features to integrate 
  id.features <- SelectIntegrationFeatures(object.list = seurats, nfeatures = nrow(seurats[[1]]))
  
  # make a reference list for the next step  
  anchors <- FindIntegrationAnchors(object.list = seurats, dims = 1:30)
  
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
  
  # now we add the cell type info from cell assigh to the metadata, we start by extracting this info from the sces 
  cell_types <- lapply(sces, function(sce) sce$cell_type)
  cell_types <- as.list(c(unlist(list(cell_types))))  # prepare for rbind by putting data in a good format
  cell_types <- do.call(rbind, cell_types)
  
  # add cell types to metadata 
  # integrated <- AddMetaData(integrated, cell_types, col.name = 'cell_types') # this one doesnt work for some reason, but code below does
  integrated$cell_types <- cell_types
  
  # the default assay is the new one 
  DefaultAssay(integrated) <- "integrated"
  
  # then some standard workflow for visualization and dim reduction 
  integrated <- ScaleData(integrated, verbose = FALSE)
  set.seed(1998)
  integrated <- RunPCA(integrated, verbose = FALSE)
  integrated <- RunTSNE(integrated, dims = 1:30, dim.embed = 3, seed.use = 300)
  integrated <- RunUMAP(integrated, dims = 1:30, seed.use = 1000)
  
  # clustering
  integrated <- FindNeighbors(integrated, dims = 1:10)
  integrated <- FindClusters(integrated, resolution = 0.5)
  
  # save the data
  saveRDS(integrated, file = integrated_object_name) # corrected data
  # saveRD(integrated, file = "/huntsman/amunzur/data/misc/important_data.gzip", compress = TRUE)
  
  # compute the combined but uncorrected object, only fix the sequencing depth 
  # load the normalized data 
  sces <- lapply(path_to_sce_cas, function(path) readRDS(path))
    
  sces <- intersect_all(sces) # subset to common genes 
  sces <- scran_batch_norm(sces)  # only correct for sequencing depth
  combined <- do.call(combine_sces, sces) # combine the objects 
    
  # some dim reduction 
  set.seed(1564)
  combined <- runPCA(combined)
  combined <- runTSNE(combined, ncomponents = 3)
  combined <- runUMAP(combined)
    
  # now save the object 
  saveRDS(combined, file = path_to_combined_sce) # uncorrected data 
    
} # end of function 

whatagreatfunction(path_to_sce_cas = args$path_to_sce_cas, 
                   path_to_combined_sce = args$path_to_combined_sce, 
                   integrated_object_name = args$integrated_object_name, 
                   ids_integration = args$ids_integration)



