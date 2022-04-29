## quick functions! ##

#this is a visualization function for seurat object to check main QC metrics i care about

preQC<-function(seur){
  seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = "^MT-")
  # Visualize QC metrics as a violin plot
  p1=VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2=FeatureScatter(seur, feature1 = "nFeature_RNA", feature2 = "percent.mt")
  plots=plot_grid(p1, p2)
  return(plots)
}


# this function to do qC on seurat based on cvisualization, allows me to dynamically change thresholds!

QC<- function(JZ,nFeature_min,nFeature_max,mt){
  JZ[["percent.mt"]] <- PercentageFeatureSet(JZ, pattern = "^MT-")
  JZ=subset(JZ, subset = nFeature_RNA > nFeature_min & nFeature_RNA < nFeature_max & percent.mt<mt)
  return(JZ)
}

# clustering and dim reuctions
NormaVar<-function(JZ, nfeatures){
  JZ<-NormalizeData(JZ)
  JZ<- FindVariableFeatures(JZ, selection.method = "vst", nfeatures =nfeatures)
  return(JZ)
}


PCA<-function(a){
  a <- ScaleData(a, features = a@assays$RNA@var.features, verbose = T)
  a <- RunPCA(a, features=a@assays$RNA@var.features, verbose = F)
  return(a)
}

Cluster<-function(a, ndims, resol){
  a <- FindNeighbors(a, dims = ndims)
  a <- FindClusters(a, resolution = resol)
  a <- RunUMAP(a, dims= ndims)
  return(a)
}


##__Enrichment Function__##
# this function takes a list of data frames of markers output from DE
run_gene_enrichment <- function(marker_list, pathways_input){
  
  pathways<-gmtPathways(pathways_input)
  # create a list of data frames for results of enrichmet per cluster/group
  fgseaRes_list <- lapply(marker_list, function(List){ 
    
    # this extracts gene symbol and logFC from every cluster
    List<-List%>%dplyr::select(gene,avg_logFC)
    List_Rank<-deframe(List)
    fgseaRes <- fgsea(pathways,List_Rank,minSize=15, maxSize=600, nperm=20000)%>%arrange(-NES)
    
    return(fgseaRes)
    
  }) 
  
  sig_pathways <- lapply(fgseaRes_list, function(x){
    
    x <- x[x$padj < 0.05 , c('pathway', 'pval', 'padj', 'NES', 'leadingEdge','size')] 
    
    return(x)
    
  })
  
  sig_pathways_name <- lapply(sig_pathways, function(x){x$pathway})
  
  output <- list(fgseaRes_list = fgseaRes_list, sig_pathways = sig_pathways, sig_pathways_name = sig_pathways_name)
  
  return(output)
  
}

##__Plotting Function__##

#gene_nerichment_list is the ouptut of the above function
make_ge_plots <- function( fgsesRes_list, pathway_name, pthresh = 0.05, ntop = 10){
  
  clusters <- names(fgsesRes_list)
  
  plots <- lapply(clusters, function(cluster, fgsesRes_list, pthresh, ntop, pathway_name, y_text_size){
    
    # select top 10 significant pathways  
    fgsesRes <- fgsesRes_list[[cluster]]
    fgsesRes <- fgsesRes[fgsesRes$padj < pthresh , ]
    fgsesRes <- top_n(fgsesRes, n = ntop, wt = abs(NES))
    
    #fgsesRes$pathway <- as.character(lapply(fgsesRes$pathway, function(pw) str_replace(pw, pattern = glue(pathway_name, "_"), "")))
    
    if(dim(fgsesRes)[1] > 0){
      
      plot <- ggplot(fgsesRes, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill=padj), width = 0.4) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title= paste("cluster", cluster, pathway_name,"pathways"))+theme_minimal()+
        theme(text = element_text(face = "bold", color = "black",size = 10),
              axis.text = element_text(face = "bold", color = "black",size =8))
      
      return(plot)
      
    } else{
      
      return(NULL)
      
    }
    
  }, fgsesRes_list= fgsesRes_list, pthresh = pthresh, ntop = ntop, pathway_name = pathway_name)
  
  plots <- plots[unlist(lapply(plots, function(plot) ! is.null(plot)))]# looks at plots list and keeps those ewntires that are not null
  # unlist
  return(plots)
  
}

# mutation_status for a seurat object:
assign_mutation<-function(seur, green, red, thresh){
  green_only <- colnames(seur)[which(as.matrix(seur@assays$RNA@data[green,] >= thresh))]
  red_only <- colnames(seur)[which(as.matrix(seur@assays$RNA@data[red,] >= thresh))]
  both <- colnames(seur)[which(as.matrix(seur@assays$RNA@data[green,] >= thresh) & as.matrix(seur@assays$RNA@data[red,] >= thresh))]
  meta<-seur@meta.data
  meta$Barcode=rownames(meta)
  cell_meta <- meta %>% dplyr::mutate(mutation_status = case_when(.$Barcode %in% setdiff(green_only, both) ~ 'green_only',
                                                                  .$Barcode %in% setdiff(red_only, both) ~ 'red_only',
                                                                  .$Barcode %in% both ~ 'both',
                                                                  TRUE ~ 'none'))
  return(cell_meta)
  }
## bar plot of comcopistions

plotCombo<-function(cell_meta,x,yy){
  cluster_meta <- cell_meta%>%dplyr::select({{x}},{{yy}})%>% group_by({{x}},{{yy}}) %>% tally(name="NumberOfCells")
  p=ggplot(cluster_meta, aes(x={{x}}, y = NumberOfCells, fill = {{yy}})) +
    geom_bar(stat = 'identity')+theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
  return(p)
}

## get normalized sce

getNormSce<-function(id){
  data <-readRDS( paste0("/huntsman/general/data/normalized/",id,"/sce_norm.rds"))
  return(data)
}

## get raw cell ranger output and create a seurat object with default qc metrics

getData<-function(path){
  data <- Read10X( paste0("/huntsman/general/data/raw/",path,"/filtered_feature_bc_matrix"))
  DH33<- CreateSeuratObject(counts =data, project = path)
  return(DH33)
}


###a function to read in the kallisto count matrix##
read_kallisto_sparse<- function(cells, regions, mtx){
  mtx<- Matrix::readMM(mtx)
  # the sparse matrix with rows are cells and columns are peaks/features
  mtx<- t(mtx)
  regions<- read_tsv(regions, col_names = FALSE)
  cells<- read_tsv(cells, col_names = FALSE)
  rownames(mtx)<- regions$X1
  # cellranger add -1 to the cell barcode, I add it for later compare with cellranger output
  colnames(mtx)<- paste0(cells$X1, "-1")
  return(mtx)
}

### this function for making a venndiagram for seurat objects that have tomato and gfp tags
ven_seurat_mutationsT<-function(seur,eGFP,dsRED){
  egfp_only <- as.vector(seur@assays$RNA@data[eGFP,] >= 0.1)
  tdTomato_only <- as.vector(seur@assays$RNA@data[dsRED,] >= 0.1)
  v <- limma::vennCounts(cbind(egfp_only, tdTomato_only))
  return(v)
}


## this function for JZ wuth egfp dsred diagram
ven_seurat_mutations<-function(seur,eGFP,dsRED){
  egfp_only <- as.vector(seur@assays$RNA@data[eGFP,] >= 0.25)
  dsred_only <- as.vector(seur@assays$RNA@data[dsRED,] >= 0.25)
  v <- limma::vennCounts(cbind(egfp_only, dsred_only))
  return(v)
}

# making boxplots from a seurat object for genes of interest 

boxplotting<-function(seur,genes, id){
  DefaultAssay(seur)="RNA"
  sub=seur[genes,]
  submeta=sub@meta.data%>%dplyr::select({{id}})
  exp=t(sub@assays$RNA@data)
  exp=data.frame(exp)
  melti=cbind(submeta,exp)
  melted=melt(melti)
  print(ggplot(melted, aes(x={{id}}, y=value, fill={{id}})) + geom_boxplot()+ theme_classic()+facet_wrap(~variable))
}



