# HOW DO I USE THIS SCRIPT? 
# all the steps in this snakefile are connected one after the other, and all the outputs and inputs are written according to 
# my file keeping system. The best thing to do would be just to update the ids and and pair_ids shown below and the run the whole analysis
# from beginning to end. These two lines are the only things you need to change: 
# However, as of March 17, note that you can only write one pair in the pair_ids part below. So if you need to run more that one pair, you 
# need to do it one pair at a time. 

# most of the scripts we have will use these ids. regardless of the analysis you do, write all the samples you include in your analysis. 
# we use the ids here for making sces, qc, normalization and clustering. 
ids = ['DH3', 'DH4', 'DH10', 'DH15', 'DH16', 'DH17', 'VOA11068_ENOC', 'DH13', 'DH18']
pair_ids = ['DH15-DH16']

# LOOK HERE FOR DGE:
# before we can do dge, we need to integrate them to remove batch effects. separate replicates by '-'. DO ONE GROUP AT A TIME.
ids_integration = ['DH3-DH4-DH10-DH15-DH16-DH17']

# if you are doing DGE analysis, write your pairs here. separate replicates by '-', separate different treatments by '='. 
# second pair will be compared to the first one. 
ids_dge = ['DH3-DH4-DH10-DH15-DH16-DH17=VOA11068_ENOC-DH13-DH18']

# once you run the pipeline, snakemake will save all the outputs in my folders. there is no need to change the following paths. 

# these are the locations of the files snakemake will need, see below. 
# it will go through the ids list and fill the {id} slot with variables from ids. 
# it will do that for each id in the ids list 

# MODIFYING SCES: 

# make sces
sce_raw = expand('../data/processed/{id}/sce.rds', id = ids)

# perform quality control
sce_qc = expand('../data/qc/{id}/sce_qc.rds', id = ids)

# normalize 
sce_norm = expand('../data/normalized/{id}/sce_norm.rds', id = ids)

# clustering 
sce_clus = expand('../data/clustered/sce/{id}/sce_clus.rds', id = ids)

# combined, but not integrated (has batch effects)
sce_uncorrected = expand('../data/integrated/{ids_integration}/uncorrected.rds', ids_integration = ids_integration)

# integration (batch effects removed)
seurat_integ = expand('../data/integrated/{ids_integration}/integrated.rds', ids_integration = ids_integration)

# MARKDOWN REPORTS, we save all of these in the reports folder (not in pipeline). 
summary_stats_report = expand('../reports/summary_stats/{pair_ids}/summary_stats.html', pair_ids = pair_ids)

# this is clustering with only 1 data set 
separate_clustering_report = expand('../reports/separate_clustering/{id}/separate_clustering.html', id = ids)

# this is the report we generate when we integrate multiple samples 
integration_report = expand('../reports/integration/{ids_integration}/integration_report.html', ids_integration = ids_integration)

# COMBINED CLUSTERING 
# where you save the output report 
combined_clustering_report = expand('../reports/combined_clustering/{pair_ids}/combined_clustering.html', pair_ids = pair_ids)
# where you save the uncorrected data 
combined_clustering_uncorrected = expand('../data/clustered/combined/{pair_ids}/uncorrected_seurat_object.rds', pair_ids = pair_ids)
# where you save the batch corrected data
combined_clustering_integrated = expand('../data/clustered/combined/{pair_ids}/integrated_seurat_object.rds', pair_ids = pair_ids)

# paired differential gene expression analysis 
paired_dge_basic_report = expand('../reports/paired_dge_basic/{pair_ids}/paired_dge_basic.html', pair_ids = pair_ids)

# DGE analysis using edgeR 
edgeR_report = expand('../reports/dge_edgeR/{ids_dge}/dge_analysis.html', ids_dge = ids_dge)

# compare scran and edgeR 
edgeR_results_up = expand('../data/dge/edgeR_up/{pair_ids}', pair_ids = pair_ids) # upregulated genes in the transduced sample edgeR found 
edgeR_results_down = expand('../data/dge/edgeR_down/{pair_ids}', pair_ids = pair_ids) # downregulated genes in the transduced sample edgeR found 
compare_dge_report = expand('../reports/dge_comparison/{pair_ids}/dge_comparison.html', pair_ids = pair_ids)   # this is where we save the knitted html file that includes the comparison data 

# integrate multiple data sets 
rule all:
  input:
    # sce_raw,
    # sce_qc,
    # sce_norm,
    # sce_clus,
    # sce_uncorrected, # only combined sces, without removing batch effects 
    # seurat_integ, # combined sces with batch effects removed 
    # integration_report, # dim reduction plots before and after batch effect removal 
    # separate_clustering_report, # cluster an sce alone
    # summary_stats_report,
    edgeR_report,
    # paired_dge_basic_report

# convert 10X counts to an sce 
rule make_sce:
  input:
    "../data/raw/{id}"
  output:
    "../data/processed/{id}/sce.rds"
  params:
    "{id}"
  shell:
    "Rscript pipeline/data_preparation/convert_to_sces.R \
     --path_to_10X {input} \
     --output_file_name {output} \
     --id {params} "


# perform quality control on the sce 
rule make_sce_qc:
  params:
    whichMethod = 'default', # to use different parameters here, just update them.
    min_features = 1000,
    mito_thresh_max = 15,
    mito_thresh_min = 1,
    nmads = 3,
    seed = 1756
  input:
    '../data/processed/{id}/sce.rds'
  output:
    '../data/qc/{id}/sce_qc.rds'
  shell:
    "Rscript pipeline/qc/default_qc.R \
     --whichMethod {params.whichMethod} \
     --path_to_sce {input} \
     --output_file_name {output} \
     --mito_thresh_max {params.mito_thresh_max} \
     --mito_thresh_min {params.mito_thresh_min} \
     --nmads {params.nmads} \
     --seed {params.seed} \
     --min_features {params.min_features}"


# normalize the sce by taking log counts and perform dim reduction at the same time 
rule normalize_sce:
  params:
    seed = 1756
  input:
    '../data/qc/{id}/sce_qc.rds'
  output:
    '../data/normalized/{id}/sce_norm.rds'
  shell:
    "Rscript pipeline/normalization/normalize_sce.R \
     --path_to_QCed_sce {input} \
     --output_file_name {output} \
     --seed {params.seed} "

# cluster the sce alone 
rule cluster_sce:
  params:
    curr_dir = os.getcwd(),
    sample = '{id}'
  input:
    '../data/normalized/{id}/sce_norm.rds'
  output:
    output_path = '../data/clustered/sce/{id}/sce_clus.rds', # this is where we save the sce with clustering data
    report = '../reports/separate_clustering/{id}/separate_clustering.html' # this is where we save the html report
  shell:
     "Rscript -e \"rmarkdown::render('pipeline/clustering/separate/clustering.Rmd',\
     output_file='{params.curr_dir}/{output.report}', \
     knit_root_dir='{params.curr_dir}',\
     params = list(id ='{params.sample}', output_path = '{output.output_path}'))\" "


# integrate (remove batch effects) and cluster 2 sces together 
rule combined_clustering:
  params:
    curr_dir = os.getcwd(), # project directory, parent of this Snakefile
    pair_ids = pair_ids,
  input:
    sce_clus
  output:
    report = combined_clustering_report,
    output_uncorrected = combined_clustering_uncorrected,
    output_integrated = combined_clustering_integrated
  shell:
    "Rscript -e \"rmarkdown::render('pipeline/clustering/combined/combined_cluster.Rmd',\
     output_file='{params.curr_dir}/{output.report}', \
     knit_root_dir='{params.curr_dir}',\
     params = list(ids ='{params.pair_ids}', output_uncorrected = '{output.output_uncorrected}', output_integrated = '{output.output_integrated}'))\" "

# structure of shell goes like this in each line:
# markdown file we are running
# where we save the output, the html file
# where the markdown file is running. just put the current directory
# these are the parameters


# make summary stats for two sces at the same time 
rule make_summary_stats:
  params:
    curr_dir = os.getcwd(),
    id = ids,
    pair_ids = pair_ids
  input:
    sce_raw,
    sce_qc,
    sce_norm
  output:
    output_path = '../data/integrated/{id_list}/integrated.rds'
  shell:
    "Rscript -e \"rmarkdown::render('pipeline/summary_stats/summary_stats.Rmd',\
     output_file='{params.curr_dir}/{output.report}', \
     knit_root_dir='{params.curr_dir}',\
     params = list(ids ='{params.pair_ids}'))\" "
  
# remove batch effects in multiple sces,      
# we use this next rule for integrating more than 2 samples. note that this rule just generates the integrated objects, it doesnt 
# make the integration report. we have the next rmarkdown rule for that. 
rule multiple_integration: 
  input:
    sce_qc, # adding these here as inputs helps snakemake understand that we need to have QCed and normalized sces before we can run this rule 
    sce_norm
  output:
    output_path_uncorrected = sce_uncorrected,
    output_path_integrated = seurat_integ
  shell:
    "Rscript pipeline/integration/integration.R \
     --path_to_sce_qc {input} \
     --path_to_sce_norm {input} \
     --output_file_name_uncorrected {output} \
     --output_file_name_integrated {output} \
     --ids_integrated {input}"
     
# show dim reduction plots before and after batch effect removal. 
rule integration_report: 
  params: 
    curr_dir = os.getcwd(),
    ids = ids_integration
  input:
    sce_uncorrected, # these objects will be computed by the multiple integration rule, just above  
    seurat_integ 
  output:
    report = integration_report
  shell:
    "Rscript -e \"rmarkdown::render('pipeline/integration/integration_report.Rmd',\
     output_file='{params.curr_dir}/{output.report}', \
     knit_root_dir='{params.curr_dir}',\
     params = list(ids ='{params.ids}'))\" "

# compute dge analysis across multiple samples using edgeR 
rule edgeR_basic:
  params:
    curr_dir = os.getcwd(),
    ids = ids_dge
  input:
    sce_uncorrected,
    seurat_integ
  output:
    report = edgeR_report
    # output_path_up =
    # output_path_down =
  shell:
    "Rscript -e \"rmarkdown::render('pipeline/paired_dge/edgeR_basic.Rmd',\
     output_file='{params.curr_dir}/{output.report}', \
     knit_root_dir='{params.curr_dir}',\
     params = list(ids ='{params.ids}'))\" "
    
# compute dge analysis in two samples using scran      
rule paired_dge_basic: 
  params: 
    curr_dir = os.getcwd(),
    ids = pair_ids
  input:
    sce_qc 
  output:
    output_path_paired = '../data/dge/paired_dge_basic/{pair_ids}/outputsTEST.rds',
    report = '../reports/paired_dge_basic/{pair_ids}/paired_dge_basic.html'
  shell:
    "Rscript -e \"rmarkdown::render('pipeline/paired_dge/paired_dge_basic.Rmd',\
     output_file='{params.curr_dir}/{output.report}', \
     knit_root_dir='{params.curr_dir}',\
     params = list(ids ='{params.ids}, {output.output_path_paired}'))\" "


