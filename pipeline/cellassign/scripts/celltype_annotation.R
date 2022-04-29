# Celltype annotation script for the DBZ paper
# 
# Purpose:
# To provide an alternative method by which to identify celltype than cellassign
# This will be done by clustering and then subsequently manually annotating clusters.
# This script is designed to work hand-in-hand with /huntsman/minhbui/DBZ/pipeline/paper/supplementary_plots.R
# May 11th, 2020

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(viridis)
  library(BiocParallel)
  library(Seurat)
})

get_sces <- function(ids){
  sces <- lapply(ids, function(id) readRDS(paste('/huntsman/general/data/clustered/sce/', id, '/sce_clus.rds', sep = '')))
  names(sces) <- ids
  return(sces)
}

#########################################################
# Create a Heatmap
#########################################################

# Enter sample IDs
#ids <- c('VOA10819UT_DBZ', 'VOA10819UT_control', 'VOA10286UT_DBZ', 'VOA10286UT_control')
#ids = c('DH1_DBZ', 'DH1_control', 'DH21', 'DH21_control', 'DH22', 'DH22_control', 'DH23', 'DH32', 'DH35')
#ids = c('VOA10286UT_control', 'DH11_control', 'DH22_control')
ids = c('VOA10819UT_control', 'DH21_control')

# Load sample data from clustered/sce. If not present, run dimred_cluster on Snakefile
sces <- get_sces(ids)

# Load marker lists
oldMarkers = read.csv("/huntsman/mdouglas/scRNAseq-workflow/pipeline/cellassign/cellassign_markers/allCellMarkers.csv", 
                      row.names = 1, col.names = c(1:50), fill = TRUE, header = FALSE, stringsAsFactors = FALSE)

markers <- read.csv("/huntsman/choran/scRNAseq-workflow/pipeline/cellassign/cellassign_markers/140322_DBZTest.csv", 
                    row.names = 1, col.names = c(1:50), fill = TRUE, header = FALSE, stringsAsFactors = FALSE)

# Original heatmap markers subset
blood <- t(cbind(oldMarkers[1,], oldMarkers[2,], oldMarkers[3,], oldMarkers[4,], oldMarkers[6,], oldMarkers[7,], oldMarkers[8,], oldMarkers[9,]))
blood <- blood[!apply(blood == "", 1, all),]
blood <- unique(trimws(unname(blood)))
rest_genes <- t(cbind(oldMarkers[5,], oldMarkers[10,]))
rest_genes <- rest_genes[!apply(rest_genes == "", 1, all),]
rest_genes <- unique(trimws(unname(rest_genes)))
features <- c(blood, rest_genes)
markerType = "OrigSubset"

# Original heatmap markers ALL
markerList <- t(cbind(oldMarkers[1,], oldMarkers[2,], oldMarkers[3,], oldMarkers[4,], oldMarkers[6,], oldMarkers[7,], oldMarkers[8,], oldMarkers[9,], oldMarkers[10,], oldMarkers[11,], oldMarkers[12,], oldMarkers[13,], oldMarkers[14,]))
markerList <- markerList[!apply(markerList == "", 1, all),]
features <- unique(trimws(unname(markerList)))
markerType = "OrigAll"

# Major cell type Markers
markerList = t(cbind(markers[33,], markers[34,], markers[35,], markers[36,]))
markerList <- markerList[!apply(markerList == "", 1, all),]
features <- unique(trimws(unname(markerList)))
markerType = "MajorType"

# Slideshow Markers
markerList = t(cbind(markers[21,], markers[22,], markers[23,], markers[24,], markers[26,]))
markerList <- markerList[!apply(markerList == "", 1, all),]
features <- unique(trimws(unname(markerList)))
markerType = "Slideshow"

# Maya Markers
markerList = t(cbind(markers[21,], markers[22,]))
markerList <- markerList[!apply(markerList == "", 1, all),]
features <- unique(trimws(unname(markerList)))
markerType = "Maya"

# If you only want to cluster by listed markers and not all genes in samples
markerOnlySCE = list()
# Remove all genes except markers and run new dimred
for(ind in 1:length(sces)) {
  geneIndex = which(rownames(assays(sces[[ind]])$counts) %in% features)
  markerOnlySCE[[ind]] = sces[[ind]][geneIndex,]
  markerOnlySCE[[ind]] <- suppressWarnings(runPCA(markerOnlySCE[[ind]]))
}

# Cluster cells by PCA
for(i in seq(1:length(sces))){
  sce <- markerOnlySCE[[i]]
  g <- buildSNNGraph(sce, k=50, use.dimred = 'PCA', BPPARAM = MulticoreParam(workers = 12))
  clust <- igraph::cluster_walktrap(g)$membership
  sce$cluster <- factor(clust)
  markerOnlySCE[[i]] <- sce
}


for(sample_num in 1:length(sces)) { 
  sce <- markerOnlySCE[[sample_num]]
  seur_sce <- as.Seurat(sce)
  Idents(seur_sce) <- seur_sce@meta.data[["cluster"]]
  seur_sce <- ScaleData(seur_sce)
  
  H3 <- DoHeatmap(seur_sce, 
                  features = features,
                  size = 3
       ) + scale_fill_viridis() + theme(legend.position="none")
  
  
  ggsave(H3, width = 12, height = 12, filename = paste('/huntsman/choran/scRNAseq-workflow/pipeline/cellassign/scripts/heatmaps/NEW', ids[sample_num], markerType, '_markers.png', sep = ''))
}


id_list = c(markerOnlySCE[[1]]@colData@listData[["id"]], markerOnlySCE[[2]]@colData@listData[["id"]])
type_list = c(markerOnlySCE[[1]]@colData@listData[["cluster"]], markerOnlySCE[[2]]@colData@listData[["cluster"]])

celltype_df = data.frame(id_list, type_list)
colnames(celltype_df) = c("id", "type")
sum(celltype_df[celltype_df$id == "DH21_control",]$type %in% c(1,2,3,4))

celltype_table = data.frame(c("EO2", "EO2", "EO1", "EO1", "EO3", "EO3", "EO4", "EO4", "EO5", "EO5"), c("Ciliated", "Unciliated", "Ciliated", "Unciliated", "Ciliated", "Unciliated", "Ciliated", "Unciliated", "Ciliated", "Unciliated"), c(366, 1476, 1274, 1956, 403, 3621, 0, 3483, 0, 2700))
colnames(celltype_table) = c("Experiment", "Type", "Proportion")

stacked_bar = ggplot(celltype_table, aes(fill=Type, y=Proportion, x=Experiment)) + 
  scale_fill_viridis(option = "plasma", discrete = T) +
  geom_bar(position="fill", stat="identity")

stacked_bar = ggplot(celltype_table, aes(fill=Type, y=Proportion, x=Experiment)) + 
  scale_fill_manual(values=c("#68116B", "#F229FA")) +
  geom_bar(position="fill", stat="identity")

stacked_bar

specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)



