# scRNAseq-workflow

> This repository contains the code used by the Huntsman Lab at BCCRC for single cell RNA sequencing analysis from 10x experiments.

**The steps of the pipeline are as follows:**  
1. process_qc_normalization: Make a SingleCellExperiment object from filtered counts. Perform quality control to further filter cells with low reads and/or high mitochondrial content. FInally, we perform normalization. We run this rule to get the normalized sce object. You can also set save_qc to "yes" if you want the qc'ed sce object.
2. run_summary_stats: Use the Seurat workflow to visualize the results of data before and after quality control. Results include summary statistics, HVGs (volcano plots), clustering and dimension reduction plots, diffusion maps, and heatmaps showing signature genes for each cluster. We run this rule to generate a rmd report on summary statistics and dimension reduction plots.
3. dimred_cluster: Do dimensionality reduction (PCA, tSNE and UMAP). Perform unsupervised clustering with scran/scater. We run this rule to get the sce object with clustering and dimension reduction information.
4. compute_diff_map: Compute diffusion map coordinates. We no longer run this rule - instead, we compute diffusion map coordinates and do the visualization in run_summary_stats rule.
5. do_batch_correction: Do batch correction and integration using seurat and scran if we have more than one sample. We run this rule to get scran corrected and seurat corrected objects for downstream analyses (i.e., integrated cell assign); and generate a rmd report visualizing dimension reduction and clustering plots before and after batch correction. 
6. integrate_cell_assign_results: Annotate cell types using cell assign and visulize cell clustering and cell type with dimension reduction plots. Although it is called integrate_cell_assign, it can handle both one-sample (with no integration) and two-sample (with integration) cases. Supply the argument `cell_type` with cell types you want to include for your sample. Choose "no" for `integration` if you don't want to perform the integration step (e.g., if running one sample). We run this rule to save a csv file on cell type information; and generate a rmd report. 
7. run_cell_cycle_report: Assign cell cycle phases using cyclone and visualize cell cycle phases with dimension reduction plots and bar plots. We run this rule to save a csv file on cell cycle phase information; and generate a rmd report. 
8. DGE_scran, DGE_edgeR_TWO_samples, DGE_edgeR_MULTIPLE_samples: Calculate the differentially expressed genes using scran (two paired samples) or edgeR (multiple samples), visualize the result with volcano plots, and perform enrichment analysis to find upregulated pathways.
9. seurat_to_loom: convert sce with cell type information to seurat and then to a loom file in preparation for velocity analyses. This loom file contains information on cell clustering and/or cell type information (if we have run cell assignment). We run thsi rule to get a loom object.

For workflow management, we use Snakemake. For details on how to install and use it, refer <a href="https://snakemake.readthedocs.io/en/stable/" target="_blank">here</a>. More detailed instructions on how to run the pipeline through Snakemake are given in the Snakefile, and in the snakefile help document as well. 

#### Design of the Snakefile  
Lists that we use as wildcards are given at the top of the Snakefile. When you have new samples to analyze, you need to add the samples you want to investigate. Details of each list will be explained below. Next, paths where the computed objects will be saved are given. You can leave them as they are, or if you would like to save your data elsewhere, you can change them as well. If you choose to change them, make sure the path you supply doesn't end with '/'. For example, supplying `'../data/clustered/sce/{samples}/sce_clus.rds/'` would be wrong; instead you should supply `'../data/clustered/sce/{samples}/sce_clus.rds'`. the last two The last thing the user needs to pay attention to is the "rule all" section of the Snakefile, as shown here:  
```
rule all:  
  input:  
    # process_qc_normalization,
    # run_summary_stats,
    # dimred_cluster,
    # do_batch_correction,
    # integrate_cell_assign_results,
    # run_cell_cycle_report,
    # seurat_to_loom,
    
    # DGE_scran,
    # DGE_edgeR_TWO_samples,
    # DGE_edgeR_MULTIPLE_samples,

    ## compute_diff_map, # we now compute and save dm coordinates in run_summary_stats step
    ## perform_dim_reduction, # determine number of PCs to use
    ## cluster_sce, #clustering on dividual sample
    ## cluster_two_samples, #clustering on integrated sample
    ## find_common_dges, # common differentially expressed genes in a paired sample
```  

User should uncomment the name of the analyses they wish to do to, there is not a limit on the number of analyses that can be done. Pipeline can be run with the `snakemake` command. This command should be typed to the terminal window in RStudio.   

**Removing mitochondrial and ribosomal genes**  
When doing quality control, you have the option to remove mitochondrial and ribosomal genes. If you wish to remove them, navigate to the "process_qc_normalization" rule and "run_summary_stats" rule in the Snakefile and set the "remove_mito_and_ribo" parameter to "yes". If you write "no", mitochondrial and ribosomal genes will be retained.  

**Subsetting the SCE object to HVG**  
In the dim step, you also have the option to subset the sces to highly variable genes (HVG). If you set the "HVG" parameter in the "dim_cluster" rule, you also need to need to indicate a percentage in the "top_HVG" parameter. The percentage you write here will be used to determine the percentage of genes with highest biological components. Note that pipeline will give an error if you write a fraction here, you must supply an integer. Default number provided is 20. SCE object will subsetted to HVG, but original count matrix will be retained as well as an alternate expression named "original". You can execute the following code to get the original counts with all the genes present:  
```
# to recover original counts 
sce_qc_original <- altExp(sce_qc_hvg, "original", withColData=TRUE)
```

### Running the pipeline through Snakemake 
The following steps of the scRNA workflow use the same wildcards:  
- make a SingleCellExperiment object; perform quality control; log normalization  
- make summary statistics report
- dimensionality reduction  
To perform these analyses, fill in the ids in the Snakefile as shown:  
`ids = ['DH4', 'DH17', 'DH10', 'DH3', 'DH15', 'DH16']`  
You can write as many as you need to; analysis steps mentioned above will be computed for each of the ids; however, some analyses may treat the ids differently. For example, quality control will be performed on each id individually, but combined clustering analysis will combine all the samples given. Note that you can only do combined clustering analysis on two samples. Individual clustering and combined clustering will generate reports, which can be found in `/reports/separate_clustering` and `/reports/combined_clustering`. Once you are done with these steps, you may move on to the next analyses which require more than one id to perfom.  

Some of the analyses mentioned above requires user to pass some parameters, mentioned below:  

#### QUALITY CONTROL (make_sce_qc)  
Performs quality control of the processed SCE based on parameters passed:  
  
**whichMethod:** pass "default" or "quantile".  
"default" will filter cells using threshold for mitochondrial and ribosomal content. To filter cell sizes based on library size and gene content, outlier cells are detected based on median-absolute-deviation (MAD).  
  
**min_features:** pass a numeric value  
Cells that have less genes detected than this number will be removed. Default value is 1000.  
  
**mito_thresh_max:**  pass a whole number between 1 and 100  
Cells with mitochonrial percentage that is a higher than this number will be removed. Default value is 25.  
  
**mito_thresh_min:** pass a whole number between 1 and 100  
Cells with mitochonrial percentage that is lower than this number will be removed. Default value is 2.  
  
**ribo_thresh_max:**  pass a whole number between 1 and 100  
Cells with ribosomal percentage that is higher than this number will be removed. Default value is 60.  
  
**nmads:** pass a whole number  
This number if used to find outlier cells based on median-absolute-deviation (MAD). Default value is 3.  
  
**seed:** pass any number  
Seed is used for reproducibility.  
  
**remove_mito_and_ribo:** pass "yes" or "no", with quotation marks.  
"yes" will remove all mitochondrial and ribosomal genes, "no" will retain them. Default is set to "yes".  

**save_qc:** pass "yes" or "no", with quotation marks.
Sometimes, we might want to obtain the qc'ed sce object. You can set this option to "yes", and the pipeline will save sce_qc.rds in `/data/qc`. You should keep this setting to "no" for the sake of saving space. 

#### DIMENSIONALITY REDUCTION (dim_reduction)  
Performs PCA, t-SNE and UMAP approximations on the SCE object.  
  
**seed1:** pass any number  
**seed2:** pass any number, different from seed1  
  
**HVG:** pass "yes" or "no", with quotation marks.  
"yes" will subset the sce to highly variable genes (HVG), "no" will retain all genes.  
  
**top_HVG:** pass a number between 1 and 100  
This number is used to calculate the top percentage of genes that contribute most to variation.  
  
**top_PCs:** pass "computed" or any number  
User can pass "computed" (with quotation marks) to let the pipeline calculate the most optimal value for number of PCs to be used in PCA calculation later on. Optimal number of PCs are calculated based on the number of PCs in the denoised PCA. Note that this number will be between 5 and 50, but the user has the option to supply a specific number instead.  

#### CLUSTERING INDIVIDUAL SCES
Clusters a given sce object and projects results in various t-SNE plots.  
  
**perplexity_list:** a list of any length, including numbers between 1 and 100  
These values are used to generate t-SNE plots with various perplexities. Default values are `perplexity_list = [5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100]`. Note that t-SNE is rather robust to changes in perplexity between 5 and 50.    
  
**chosen_perp:** one single number, taken from the perplexity list given above  
After viewing t-SNE plots with various perplexities, user is expected to pick a perplexity value and save t-SNE coordinates for that value for more downstream analysis. After running the clustering script once and inspecting the t-SNE plots, you can pick one of the perplexities and run the clustering script again. The default value is set to 30.  
  
**k_list:** a list of any length, including only positive integers  
This list will be used during clustering to vary number of nearest neighbors.  You can supply any numbers, but the suggested values range between 5 and 30. Default values supplied in the pipeline are `k_list = [5, 9, 10, 12, 15, 18]`.  
  
**chosen_k:** one single number, taken from the k_list given above  
After viewing the clustering plots with various k values, the user is expected to choose a k value to continue downstream analysis. After inspecting the report once, you can pick a certain k value to store in the sce object for further downstream analysis. Default value is 15.  

**cell_types:** this is where we specify the type of cells you want to include for cellassign. The pipeline will assemble a marker gene matrix corresponding to cell types you specify. For a list of gene markers, refer to ./cellassign/scripts/create_marker_matrix.R. A spreadsheet with more details on these markers can be found here: https://docs.google.com/spreadsheets/d/1amSrlUxSD_N7XNDlqT3qvu3KLUyUYWTFux3lKH4Hnzw/edit#gid=0

  
> To run the pipeline using snakemake, change the name of the snakefile from 'sample_snakefile' to 'Snakefile' on your local machine after pulling from master.  

### Analyses including multiple data sets  

 - **Summary statistics**  
Show summary statistics for two datasets in a report. To run this, write the ids you wish to compare in the pair_ids list in the Snakefile. Make sure ids correspond to pair_ids. For example, if you are running DH22_NEW and DH22_control, you need to specify `ids=['DH22_NEW','DH22_control']`, `pair_ids=['DH22_NEW-DH22_control']`. You can supply more than one pair of ids; summary stats will be computed for each pair separately.   

 - **Batch normalization**  
Remove batch effects in two or more datasets, save the corrected samples and show findings in a report.  
To run this analysis, write down the sample names to ids_integration list shown at the beginning of the Snakefile. Separate the samples with a dash (-). Do one group of ids at a time. You also have the option to write a short explanatory message to be shown in the title of the integration report. This was designed to help readers to get a quick idea about what the report is about. You have the option to leave it blank. Write your message as a string to 'id_type' list in the Snakefile. These batch corrected data hasn't been implemented for downstream analysis.  

We use two integration methods: combine and remove batch effect using scran; integrate with seurat.
We use Seurat's integration techniques to remove batch effects between replicates. It's possible that integration causes subpopulation diversity to be eliminated as well. Run the pipeline with caution. If any of your samples have less than 200 cells, Seurat integration techniques won't work. In that case, you can try using scran batch correction, explained in more detail in the next step. Indicate ids you wish to integrate in `ids_integration` list. 

 - **DE analysis with edgeR**  
Perform DE analysis in two groups, groups can have as many replicates as you wish. Write the groups you wish to compare to 'ids_dge' list in the Snakefile. Separate replicates that belong to the same group with a '-'; separate different groups with '='. Depending on the number of samples you have, you will need to follow different paths. If you have only two or three  samples, edgeR needs to run on two samples only. You can do that by uncommenting `DGE_edgeR_TWO_samples` from the `rule all` section. If you have at least two replicates for each treatment (4 samples at least in total), uncomment `DGE_edgeR_MULTIPLE_samples` from `rule all`. Remember that second if you supply wilk compared to the first id. For example, if you supply `ids_dge = [DH1-DH2=DH3-DH4]`, DH3 and DH4 will be compared to DH1 and DH2.     

 - **DE analysis with scran**  
Perform DE analysis between two ids using t tests. Write the ids to the 'pair_ids' list in the Snakefile and separate the samples with a '-'. First sample should be transduced, and second sample should be the control. Note that this DE analysis will compare the first sample to the second sample. If you have more than two samples for each to compare, it is recommended that you use the edgeR DGE method, explained below.  

 - **Cell type assignment**  
For information about how to install cell assign, refer to <a href="https://shahlab.ca/projects/cellassign/" target="_blank">shahlab.ca/projects/cellassign</a>. We recommend that you use a virtual environment to install Tensorflow. Once the installation is complete, prior to using cell assign, activate the virtual environment by running `source ~/venv/bin/activate` in the Rstudio terminal. Snakemake rule for running cell assign doesn't use any wildcards other than ids, so writing the cell ids you wish to annotate is sufficient. Cell assign run will output a few files that can be used in downstream analysis. Note that cell types will be computed individually for all the ids you supply.  

 - **Integrating and visualizing cell assign results**  
At this step, a special rule other than Seurat integration has been designed. After computing cell types, write the samples you wish to integrate to the integration_ids list at the top of the Snakefile and separate distinct samples with "-", as shown:  
`ids_integration = ['DH4-DH17-DH10']`  
You can integrate as many samples as you wish. We previously integrated/combined samples with scran and seurat in the batch correction step. Sce objects from both methods are loaded and visualized. Cell assign reports provide visualization on cell clustering and cell type clustering on individual samples and integrated samples with dimention reduction plots. Reports for individual samples are in `/reports/cellassign` (this directory contains reports from old runs as we no longer have a separate rule for running individual sample), and integrated sample reports can be found in `/reports/cellassign_integrated` (when you run integrate_cell_assign_results, output reports are in this directory - for both individual sample and multiple samples).  

Moreover, each of these steps have two types of files: one off scripts and scripts that have been integrated to the pipeline. One off scripts end with '_O.Rmd'. If you need to repeat a specific step isolated from the Snakemake workflow, you can use the one off scripts where you can specify the ids you wish to analyze.  

Packages you need to run the pipeline are listed in the sourceFiles/utilities.R folder.  

### A few Snakemake tips  
`snakemake -n` will execute a dry run.  
`snakemake -j 5` will allocate 5 cores.  
`snakemake -F` will force already computed objects to be computed again. 
