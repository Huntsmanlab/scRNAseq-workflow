

# INPUT FUNCTIONS 

# SEURAT INTEGRATION 
def seurat_integration_output_files(wildcards): 
    
    combined_sample_names = "-".join(wildcards)
    path_to_uncorrected = '../data/integrated/' + combined_sample_names + '/uncorrected.rds'
    path_to_integrated = '../data/integrated/' + combined_sample_names + '/integrated.rds'
    
    return [path_to_uncorrected, path_to_integrated]

def seurat_integration_input_files(wildcards): 
    
    i = 0
    sce_qc_list = []
    sce_norm_list = []
    
    while i < len(wildcards): 
        
        path_to_sce_qc = '../data/qc/' + wildcards[i] + '/sce_qc.rds'
        sce_qc_list.append(path_to_sce_qc)
                
        path_to_sce_norm = '../data/normalized/' + wildcards[i] + '/sce_norm.rds'
        sce_norm_list.append(path_to_sce_norm)
        
        i += 1
        
    return sce_qc_list, sce_norm_list
    
def seurat_integration_report_path(wildcards): 
  
    combined_sample_names = "-".join(wildcards)
    path_to_uncorrected = '../reports/integration/' + combined_sample_names + '/integration_report.html'
    
    return path_to_uncorrected

# COMBINED CLUSTERING REPORT 
def combined_clustering_report_input_files(wildcards):
    
    samples = wildcards[0].split('-')
    
    i = 0 
    path_list = []
    
    while i < len(samples): 
        
        path_to_sce_clus = '../data/clustered/sce/' + samples[i] + '/sce_clus.rds/'
        path_list.append(path_to_sce_clus)
        
        i += 1
        
    return path_list
    
# BATCH CORRECTION REPORT 
# wildcards are given in a list, like ['DH7', 'DH8']
def batch_correction_output_files(wildcards): 
  
  combined_sample_names = "-".join(wildcards)
  
  path_to_sce_with_batchEffects = '../data/batch_corrected/' + combined_sample_names + '/uncorrected.rds'
  path_to_sce_without_batchEffects = '../data/batch_corrected/' + combined_sample_names + '/corrected.rds'
  path_to_report = '../reports/batch_normalization/' + combined_sample_names + '/batch_correction.html'

  return [path_to_sce_with_batchEffects, path_to_sce_without_batchEffects, path_to_report]

# wildcards are given in a list, like ['DH7', 'DH8']
def batch_correction_input_files(wildcards): 
  
  i = 0
  path_to_sce_clus_list = []
  
  while i < len(wildcards): 
    
    sce_clus = '../data/clustered/sce/' + wildcards[i] + '/sce_clus.rds'
    path_to_sce_clus_list.append(sce_clus)
    
    i += 1

  return path_to_sce_clus_list


# SUMMARY STATS REPORT 

# wildcards are given like this IN A LIST: ['DH7-DH8', 'DH10-DH24']. one pair at a time 
def summary_stats_files(wildcards): 
  
  wildcards = wildcards[0].split('-')
  
  i = 0 
  paths_to_sce_raws = []
  paths_to_sce_qcs = []
  
  # paths to sces 
  while i < len(wildcards): 
    
    sce_raw = '../data/processed/' + wildcards[i] + '/sce.rds'
    sce_qc = '../data/qc/' + wildcards[i] + '/sce_qc.rds'
    
    paths_to_sce_raws.append(sce_raw)
    paths_to_sce_qcs.append(sce_qc)
    
    i += 1 
    
  return paths_to_sce_raws, paths_to_sce_qcs

# EDGER 

def edgeR_input_files(wildcards): 
    
    # separate into two groups 
    groups = wildcards[0].split('=')
    
    # further separate groups 
    sample1 = groups[0]
    sample2 = groups[1]
    
    sample1_path = '../data/integrated/' + sample1 + '/integrated.rds'
    sample2_path = '../data/integrated/' + sample2 + '/integrated.rds'
    
    paths = [sample1_path, sample2_path]
    
    # load sce clus, separate groups into individual samples
    sample1_individual = groups[0].split('-')
    sample2_individual = groups[1].split('-')
    
    # sample 1
    i = 0 
    sample1_sceclus_paths = []
    
    while i < len(sample1_individual): 
        
        sce = sample1_individual[i]
        path = '../data/clustered/sce/' + sce + '/sce_clus.rds'
        sample1_sceclus_paths.append(path)
        
        i += 1 
        
    # sample 2
    i = 0 
    sample2_sceclus_paths = []
    
    while i < len(sample2_individual): 
        
        sce = sample2_individual[i]
        path = '../data/clustered/sce/' + sce + '/sce_clus.rds'
        sample2_sceclus_paths.append(path)
        
        i += 1 
        
    sce_clus_paths = [sample1_sceclus_paths, sample2_sceclus_paths]
    
    return [paths, sce_clus_paths]
        
# PAIRED DGE BASIC 
# wildcards are given in a list like this: ['DH1-DH2', 'DH3-DH4']
def paired_dge_input_files(wildcards): 
    
    separate_samples = wildcards[0].split('-')
    inputFile_sample1 = '../data/clustered/sce/' + separate_samples[0] + '/sce_clus.rds'
    inputFile_sample2 = '../data/clustered/sce/' + separate_samples[1] + '/sce_clus.rds'
    
    return [inputFile_sample1, inputFile_sample2]


# INTEGRATE CELL ASSIGN 
# wildcards are given as a string in a list, separated by '-', like this: ['DH1-DH2']
def integrate_cell_assign_files(wildcards): 
  
  separate_samples = wildcards[0].split('-')
  
  i = 0 
  sce_cas_path_list = [] # input files 
  
  while i < len(separate_samples): 
    
    sce_cas_path = '../data/cellassign/' + separate_samples[i] + '/sce_norm_cas.rds'
    sce_cas_path_list.append(sce_cas_path)
    
    i += 1 
    
  path_to_integrated_output = '../data/integrated_cellassign/' + wildcards[0] + '/integrated_cas.rds'
  
  return sce_cas_path_list, path_to_integrated_output

# generate report for the integrated sample

def cell_assign_multiple_report_files(wildcards): 
    
    separate_samples = wildcards[0].split('-')
        
    i = 0 
    sce_cas_file_paths = []
    
    while i < len(separate_samples): 
        
        sce_cas_path = '../data/cellassign/' + separate_samples[i] + '/sce_norm_cas.rds'
        sce_cas_file_paths.append(sce_cas_path)
        
        i += 1 
        
    path_to_uncorrected = '../data/integrated/' + wildcards[0] + '/uncorrected.rds'
    path_to_integrated = '../data/integrated_cellassign/' + wildcards[0] + '/integrated_cas.rds'
    
    return sce_cas_file_paths, path_to_uncorrected, path_to_integrated
  
  
  
  
  













