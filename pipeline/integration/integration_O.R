# This script is for integrating multiple data sets. integrating allows us to combine data sets and remove batch effects at the same time 
# we may need to integrate data in various situations, so this script can be used in multiple types of analysis, such as: 
# running a clustering analysis on the combined sample 
# run dim reduction plots on combined plots 
# run DGE analysis across samples, we will need to remove batch effects for that. we can integrate data sets first here, then split if needed. 
# data here is saved into the data/integrated folder/<name of integrated samples>
# this is the one off script, we also have a version integrated into snakemake 

# load a few necessary things 
library(here)
source(here('..', 'sourceFiles', 'utilities.R'))

# update here when you have a new list. update the same update the 
id.list <-list('DH4', 'DH17')
id.orig <- id.list

# lets start! 
# load the data, note that we aren't using the normalized data since we will use seurat's normalization method 
sces <- lapply(id.list, function(id) readRDS(here('..', 'data', 'qc', id, 'sce_qc.rds')))

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

# the default assay is the new one 
DefaultAssay(integrated) <- "integrated"

# then some standard workflow for visualization and dim reduction 
integrated <- ScaleData(integrated, verbose = FALSE)
set.seed(1998)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- RunTSNE(integrated, dims = 1:30)


# save the data, assuming you made the folder where the data would be saved 
saveRDS(integrated, file = here('..', 'data', 'integrated', paste(unlist(id.orig), collapse = '-'), 'integrated.rds' )) # corrected data 

# UNCORRECTED DATA
# here we do the combined but uncorrected sample, with the batch effects
# load the normalized data 
sces_norm <- lapply(id.orig, function(id) readRDS(here('..', 'data', 'normalized', id, 'sce_norm.rds')))
sces_norm <- intersect_all(sces_norm) # subset to common genes 

# combine the normalized sces 
uncorrected <- do.call(combine_sces, sces_norm)

# some dim reduction 
set.seed(1564)
uncorrected <- runPCA(uncorrected)
uncorrected <- runTSNE(uncorrected)
uncorrected <- runUMAP(uncorrected)

# now save the object 
saveRDS(uncorrected, file = here('..', 'data', 'integrated', paste(unlist(id.orig), collapse = '-'), 'uncorrected.rds' )) # uncorrected data 


















