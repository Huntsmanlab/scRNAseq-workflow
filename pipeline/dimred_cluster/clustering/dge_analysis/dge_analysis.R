
# source a function from scran 
source('/huntsman/general/scRNAseq-workflow/pipeline/sourceFiles/getMarkerEffects_scran.R')

# load the data 
readRDS()

# here we do pairwise comparisons between clusters for each gene, and return a list of dfs containing ranked candidate markers for each cluster.
# for example, 8th element in this list would give you the candidate genes for cluster 8. 
markers.sce <- findMarkers(sce_clus, sce_clus$cluster)

# pick a cluster to further investigate 
chosen <- "1"

# subsets the interesting markerd
interesting <- markers.sce[[chosen]]
colnames(interesting)

interesting[1:10,1:4]
best.set <- interesting[interesting$Top <= 6,]
logFCs <- getMarkerEffects(best.set)

pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))

# find upregulated genes only: 
markers.sce.up <- findMarkers(sce, sce$cluster, direction="up")
interesting.up <- markers.sce.up[[chosen]]
interesting.up[1:10,1:4]




markers.sce.up2 <- findMarkers(sce, sce$cluster, 
                                direction="up", lfc=1)
interesting.up2 <- markers.sce.up2[[chosen]]
interesting.up2[1:10,1:4]

best.set <- interesting.up2[interesting.up2$Top <= 5,]
logFCs <- getMarkerEffects(best.set)

pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))


