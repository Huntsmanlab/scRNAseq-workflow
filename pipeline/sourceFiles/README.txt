This file has information about the functions used in the utilities file. 

1. visualize_cellassign()
As an input, it takes a seurat object with cell_type information and generates a dim reduction plot showing different cell types. Function will always give the same colors for the cell types.

 - seurat_object: your seurat object with cell type information 
 - reduction_type: 'tsne' or 'umap'. Must be all lower case. Must be a string.
 - group_by: which groups to highlight, either 'cell_types' or 'id'
 - master_color_palette: your color palette to choose color from 
 - master_cell_types: all the cell types cell assign can detect

master_color_palette and master_cell_types are given in the utilities file. 

2. intersect_all()
Given a list of sces, find common genes and subset all sces to the common genes. Useful during DGE analysis. 

 - sces: a list of sces, can be as long as you like
 
3. seuratClustering()
Given an sce, turn it into a seurat object and perform seurat style clustering

 - sce: an sce
 - id: the sample name of the sce to add to the seurat object 
 
4. combine_sces2()
This is a useless function, but I have to keep it here in case one of my scripts call this silly function. 
Given a list of seurat objects, combine them. 

 - seurat_list: a list of seurat objects you wish to merge 
 - ids_list: the ids of those seura objects 
 
5. make_upreg_table()

USED WITH EDGER

Given qlf (quasi-likelihood F test, an output of edgeR), make a table of upregulated genes. 
Genes are ordered in decreasing log fold change (highest to lowest) and increasing p values (lowest to highest). The final column named 'summed ranks' is a way to visualize this. Genes with low summed ranks have low p values and high log fold changes. 

 - qlf: one of the outputs of edgeR 
 
6. make_downreg_table()

USED WITH EDGER

Same as above, but now this time we make a table of downregulated genes. 

7. make_pseudo_counts()
Given a list of sces, sum the expression of cells in each cluster to make pseudo bulk counts. This allows us to overcome the issue of lack of replication by treating each cluster as a sample. This function and others similar to it is commonly used in DGE analysis with edgeR because edgeR requires replication. Output given is a data table where each row is a gene and each column is a cluster.

Warning: This is a badly written function. Next few ones are better suited for pseudo-bulking. 

 - sce_list: a list of sces that you want to use in the DGE analysis
 
8. make_pseudo_counts_seurat()
Given a combined seurat object, make pseudo bulk counts for each separate sample found in the combined seurat object. This function is commonly used in various edgeR scripts. Let's say you have 4 seurat objects and you want to do DE analysis between them using edgeR. You would merge them into one seurat object and and use the merged object as an input. The output of the function is a data frame that can directly be used in edgeR to make a DGElist object. Each row is a gene, and each column is a sample. Note that when you merge the objects, you subset all the seurat objects to the common genes between the objects. 

 - sobject: a merged seurat object that contains at least two different samples
 
9. find_and_bind()
Given TWO data frames, this function will find the common rows and subset both data frames to common rows. The output is a combined data frame where the rows are the common genes and columns are the samples.   



























