# scRNAseq-workflow

> This repository contains the code used by the Huntsman Lab at BCCRC for single cell RNA sequencing analysis.

**The steps of the pipeline are as follows:**  
- make a SingleCellExperiment object from raw counts data  
- perform quality control to filter cells with low reads and/or high mitochondrial content  
- log normalization  
- dimensionality reduction (tSNE and UMAP) 
- unsupervised clustering  
- batch correction and integration, if needed  
- DGE analysis  

For workflow management, we use Snakemake. For details on how to install and use it, refer <a href="https://snakemake.readthedocs.io/en/stable/" target="_blank">here</a>. More detailed instructions on how to run the pipeline through Snakemake are given in the Snakefile, and in the snakefile help document as well. 

#### Running the pipeline through Snakemake 
The following steps of the scRNA workflow use the same wildcards:  
- make a SingleCellExperiment object  
- perform quality control  
- log normalization  
- dimensionality reduction  
- clustering  
To perform these analyses, fill in the ids in the Snakefile as shown:  
ids = ['DH4', 'DH17', 'DH10', 'DH3', 'DH15', 'DH16']  
You can write as many as you need to; analysis steps mentioned above will be computed for each of the id. Afterwards, you may perform these analyses which require more than one id to perfom:  
**seurat integration **  
Integrate two or more data sets and show findings in a report. At the momennt, we use Seurat's integration techniques to remove batch effects between replicates. It's possible that integration causes subpopulation diversity to be eliminated as well. Run the pipeline with caution. If any of your samples have less than 200 cells, this will not work. 
**batch correction**  
Remove batch effects in two or more datasets, save the corrected samples and show findings in a report. To run this analysis, write down the sample names to ids_integration list shown at the beginning of the Snakefile. Separate the samples with a dash (-). Do one group of ids at a time. You also have the option to write a short explanatory message to be shown in the title of the integration report. This was designed to help readers to get a quick idea about what the report is about. You have the option to leave it blank. Write your message as a string to 'id_type' list in the Snakefile. These batch corrected data hasn't been implemented for downstream analysis. 
**summary statistics**  
Show summary statistics for two datasets in a report. To run this, write the ids you wish to compare in the pair_ids list in the Snakefile.  
**DE analysis with edgeR**  
Perform DE analysis in two groups, groups can have as many replicates as you wish. Write the groups you wish to compare to 'ids_dge' list in the Snakefile. Separate replicates that belong to the same group with a '-'; separate different groups with '='. 
**DE analysis with scran**  
Perform DE analysis between two ids using t tests. Write the ids to the 'pair_ids' list in the Snakefile and separate the samples with a '-'. If you have more than two samples for to compare, it is recommended that you use the edgeR DGE method. 

Moreover, each of these steps have two types of files: one off scripts and scripts that have been integrated to snakemake. One off scripts end with '_O.Rmd'. If you need to repeat a specific step isolated from the Snakemake workflow, you can use the one off scripts where you can specify the ids you wish to analyze.  

















