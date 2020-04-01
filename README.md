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

For workflow management, we use Snakemake. For details on how to install and use it, refer <a href="https://snakemake.readthedocs.io/en/stable/" target="_blank">here</a>. More detailed instructions on how to run the pipeline through Snakemake are given in the Snakefile. 