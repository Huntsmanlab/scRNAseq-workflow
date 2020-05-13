# scRNAseq-workflow

> This repository contains the code used by the Huntsman Lab at BCCRC for single cell RNA sequencing analysis from 10x experiments.

**The steps of the pipeline are as follows:**  
1. Make a SingleCellExperiment object from filtered counts  
2. Perform quality control to further filter cells with low reads and/or high mitochondrial content  
3. Visualize the results of quality control in summary stats 
4. Log normalize
5. Do dimensionality reduction (tSNE and UMAP) 
6. Perform unsupervised clustering  
7. Do batch correction and integration, if needed  
8. Calculate the differentially expressed genes and visualize the analysis  
9. Annotate cell types using cell assign 

For workflow management, we use Snakemake. For details on how to install and use it, refer <a href="https://snakemake.readthedocs.io/en/stable/" target="_blank">here</a>. More detailed instructions on how to run the pipeline through Snakemake are given in the Snakefile, and in the snakefile help document as well. 

#### Running the pipeline through Snakemake 
The following steps of the scRNA workflow use the same wildcards:  
- make a SingleCellExperiment object  
- perform quality control  
- log normalization  
- dimensionality reduction  
- individual clustering (cluster one sample)  
- combined clustering (cluster multiple samples)  
To perform these analyses, fill in the ids in the Snakefile as shown:  
ids = ['DH4', 'DH17', 'DH10', 'DH3', 'DH15', 'DH16']  
You can write as many as you need to; analysis steps mentioned above will be computed for each of the ids; however, some analyses may treat the ids differently. For example, quality control will be performed on each id individually, but combined clustering analysis will combine all the samples given. Note that you can only do combined clustering analysis on two samples. Individual clustering and combined clustering will generate reports, which can be found in `/reports/separate_clustering` and `/reports/combined_clustering`. Once you are done with these steps, you may move on to the next analyses which require more than one id to perfom.  

To run the pipeline using snakemake, change the name of the snakefile from 'sample_snakefile' to 'Snakefile' on your local machine after pulling from master. 

#### Analyses including multiple data sets  
 - **Seurat integration**  
Integrate two or more data sets and show findings in a report.  
At the moment, we use Seurat's integration techniques to remove batch effects between replicates. It's possible that integration causes subpopulation diversity to be eliminated as well. Run the pipeline with caution. If any of your samples have less than 200 cells, Seurat integration techniques won't work. In that case, you can try using batch correction, explained in more detail in the next step. 
 - **Batch correction**  
Remove batch effects in two or more datasets, save the corrected samples and show findings in a report.  
To run this analysis, write down the sample names to ids_integration list shown at the beginning of the Snakefile. Separate the samples with a dash (-). Do one group of ids at a time. You also have the option to write a short explanatory message to be shown in the title of the integration report. This was designed to help readers to get a quick idea about what the report is about. You have the option to leave it blank. Write your message as a string to 'id_type' list in the Snakefile. These batch corrected data hasn't been implemented for downstream analysis.  
 - **Summary statistics**  
Show summary statistics for two datasets in a report. To run this, write the ids you wish to compare in the pair_ids list in the Snakefile.  
 - **DE analysis with edgeR**  
Perform DE analysis in two groups, groups can have as many replicates as you wish. Write the groups you wish to compare to 'ids_dge' list in the Snakefile. Separate replicates that belong to the same group with a '-'; separate different groups with '='.  
 - **DE analysis with scran**  
Perform DE analysis between two ids using t tests. Write the ids to the 'pair_ids' list in the Snakefile and separate the samples with a '-'. If you have more than two samples for to compare, it is recommended that you use the edgeR DGE method.  
 - **Cell type assignment**  
For information about how to install cell assign, refer to <a href="https://shahlab.ca/projects/cellassign/" target="_blank">shahlab.ca/projects/cellassign</a>. Once the installation is complete, prior to using cell assign, activate the virtual environment by running `source ~/venv/bin/activate` in the terminal. Snakemake rule for running cell assign doesnt use any wildcards other than ids, so writing the cell ids you wish to annotate is sufficient. Cell assign run will output a few files that can be used in downstream analysis.  
 - **Visualizing cell assign results**  
There a few options for visualizing cell assign outputs. You can either choose to visualize one sample or multiple samples together. For multiple sample visualization, samples should be integrated first. 

Moreover, each of these steps have two types of files: one off scripts and scripts that have been integrated to snakemake. One off scripts end with '_O.Rmd'. If you need to repeat a specific step isolated from the Snakemake workflow, you can use the one off scripts where you can specify the ids you wish to analyze.  

Packages you need to run the pipeline are listed in the sourceFiles/utilities.R folder. 


<a href="https://www.thesitewizard.com/" target="_blank">thesitewizard.com</a>

















