

# INPUT FUNCTIONS 

# helper function 
def remove_last_two_dirs(path): 
    path = path.split('/') # split 
    path = path[:len(path)-2] # remove last two elements 
    
    strng = '/'
    
    path = strng.join(path) + '/'
    
    return path


# SEURAT INTEGRATION 
def seurat_integration_output_files(wildcards, path_to_uncorrected, path_to_integrated): 
    
    combined_sample_names = wildcards[0]
    path_to_uncorrected = remove_last_two_dirs(path_to_uncorrected) + combined_sample_names + '/uncorrected.rds'
    path_to_integrated = remove_last_two_dirs(path_to_integrated) + combined_sample_names + '/integrated.rds'
    
    return [path_to_uncorrected, path_to_integrated]

def seurat_integration_input_files(wildcards, path_to_sceqc, path_to_scenorm): 
  
    # separate the given wildcards into separate samples 
    wildcards = wildcards[0].split('-')
    
    i = 0
    sce_qc_list = []
    sce_norm_list = []
    
    while i < len(wildcards): 
        
        path_to_sce_qc = remove_last_two_dirs(path_to_sceqc) + wildcards[i] + '/sce_qc.rds'
        sce_qc_list.append(path_to_sce_qc)
                
        path_to_sce_norm = remove_last_two_dirs(path_to_scenorm) + wildcards[i] + '/sce_norm.rds'
        sce_norm_list.append(path_to_sce_norm)
        
        i += 1
        
    return sce_qc_list, sce_norm_list
    
# def seurat_integration_report_path(wildcards): 
#   
#     combined_sample_names = wildcards[0]
#     path_to_report = '../reports/integration/' + combined_sample_names + '/integration_report.html'
#     
#     return path_to_report

# COMBINED CLUSTERING REPORT 
def combined_clustering_report_input_files(wildcards, path_to_sce_clus):
    
    path = remove_last_two_dirs(path_to_sce_clus)
    
    samples = wildcards[0].split('-')
    
    i = 0 
    path_list = []
    
    while i < len(samples): 
        
        path_to_sce_clus = path + samples[i] + '/sce_clus.rds/'
        path_list.append(path_to_sce_clus)
        
        i += 1
        
    return path_list    
    
# BATCH CORRECTION REPORT 
# wildcards are given as a string, like DH1-DH2
def batch_correction_output_files(wildcards, path_to_sce_with_batchEffects, path_to_sce_without_batchEffects, path_to_report): 
  
  combined_sample_names = "-".join(wildcards)
    
  path_to_sce_with_batchEffects = remove_last_two_dirs(path_to_sce_with_batchEffects)
  path_to_sce_without_batchEffects = remove_last_two_dirs(path_to_sce_without_batchEffects)
  path_to_report = remove_last_two_dirs(path_to_report)
  
  path_to_sce_with_batchEffects = path_to_sce_with_batchEffects + combined_sample_names + '/uncorrected.rds'
  path_to_sce_without_batchEffects = path_to_sce_without_batchEffects + combined_sample_names + '/corrected.rds'
  path_to_report = path_to_report + combined_sample_names + '/batch_correction.html'

  return [path_to_sce_with_batchEffects, path_to_sce_without_batchEffects, path_to_report]

# wildcards are given in a list, like DH1-DH2
def batch_correction_input_files(wildcards, path_to_sce_clus): 
    
  path = remove_last_two_dirs(path_to_sce_clus)
  
  wildcards = wildcards[0].split('-')
  i = 0
  path_to_sce_clus_list = []
  
  while i < len(wildcards): 
    
    sce_clus = path + wildcards[i] + '/sce_clus.rds'
    path_to_sce_clus_list.append(sce_clus)
    
    i += 1

  return path_to_sce_clus_list

# SUMMARY STATS REPORT 

# wildcards are given like this IN A LIST: ['DH7-DH8', 'DH10-DH24']. one pair at a time 
def summary_stats_files(wildcards, paths_to_sce_raws, paths_to_sce_qcs): 
  
  wildcards = wildcards[0].split('-')
    
  sce_raw_path = remove_last_two_dirs(paths_to_sce_raws)
  sce_qc_path = remove_last_two_dirs(paths_to_sce_qcs)
  
  i = 0 
  paths_to_sce_raws = []
  paths_to_sce_qcs = []
  
  # paths to sces 
  while i < len(wildcards): 
    
    sce_raw = sce_raw_path + wildcards[i] + '/sce.rds'
    sce_qc = sce_qc_path + wildcards[i] + '/sce_qc.rds'
    
    paths_to_sce_raws.append(sce_raw)
    paths_to_sce_qcs.append(sce_qc)
    
    i += 1 
    
  return paths_to_sce_raws, paths_to_sce_qcs
# EDGER 

def edgeR_input_files(wildcards, sce_integrated, sce_clus): 
    
    # separate into two groups 
    groups = wildcards[0].split('=')
    
    # further separate groups 
    sample1 = groups[0]
    sample2 = groups[1]
    
    sce_integrated_path = remove_last_two_dirs(sce_integrated)
    sce_clus = remove_last_two_dirs(sce_clus)
    
    sample1_path = sce_integrated_path + sample1 + '/integrated.rds'
    sample2_path = sce_integrated_path + sample2 + '/integrated.rds'
    
    paths = [sample1_path, sample2_path]
    
    # load sce clus, separate groups into individual samples
    sample1_individual = groups[0].split('-')
    sample2_individual = groups[1].split('-')
    
    # sample 1
    i = 0 
    sample1_sceclus_paths = []
    
    while i < len(sample1_individual): 
        
        sce = sample1_individual[i]
        path = sce_clus + sce + '/sce_clus.rds'
        sample1_sceclus_paths.append(path)
        
        i += 1 
        
    # sample 2
    i = 0 
    sample2_sceclus_paths = []
    
    while i < len(sample2_individual): 
        
        sce = sample2_individual[i]
        path = sce_clus + sce + '/sce_clus.rds'
        sample2_sceclus_paths.append(path)
        
        i += 1 
        
    sce_clus_paths = [sample1_sceclus_paths, sample2_sceclus_paths]
    
    type(paths)
    type(sce_clus_paths)
    
    return [paths, sce_clus_paths]
    
    
def edgeR_input_files_TWO_samples(wildcards, sce_clus): 
  
    path_to_sce_clus = remove_last_two_dirs(sce_clus)
   
    # separate into two groups 
    groups = wildcards[0].split('=')
    
    # further separate groups 
    sample1 = groups[0]
    sample2 = groups[1]
    
    # load sce clus, separate groups into individual samples
    sample1_individual = groups[0].split('-')
    sample2_individual = groups[1].split('-')
    
    # sample 1
    i = 0 
    sample1_sceclus_paths = []
    
    while i < len(sample1_individual): 
        
        sce = sample1_individual[i]
        path = path_to_sce_clus + sce + '/sce_clus.rds'
        sample1_sceclus_paths.append(path)
        
        i += 1 
        
    # sample 2
    i = 0 
    sample2_sceclus_paths = []
    
    while i < len(sample2_individual): 
        
        sce = sample2_individual[i]
        path = path_to_sce_clus + sce + '/sce_clus.rds'
        sample2_sceclus_paths.append(path)
        
        i += 1 
        
    return sample1_sceclus_paths, sample2_sceclus_paths
    
        
# PAIRED DGE BASIC 
# wildcards are given in a list like this: ['DH1-DH2', 'DH3-DH4']
def paired_dge_input_files(wildcards, sce_clus): 
    
    path = remove_last_two_dirs(sce_clus)
    
    separate_samples = wildcards[0].split('-')
    inputFile_sample1 = path + separate_samples[0] + '/sce_clus.rds'
    inputFile_sample2 = path + separate_samples[1] + '/sce_clus.rds'
    
    return [inputFile_sample1, inputFile_sample2]


# INTEGRATE CELL ASSIGN 
# wildcards are given as a string in a list, separated by '-', like this: ['DH1-DH2']
def integrate_cell_assign_files(wildcards, sce_norm_cas, integrated_cas ): 
    
  path_to_cas = remove_last_two_dirs(sce_norm_cas)
  path_to_integ = remove_last_two_dirs(integrated_cas)
  
  separate_samples = wildcards[0].split('-')
  
  i = 0 
  sce_cas_path_list = [] # input files 
  
  while i < len(separate_samples): 
    
    sce_cas_path = path_to_cas + separate_samples[i] + '/sce_norm_cas.rds'
    sce_cas_path_list.append(sce_cas_path)
    
    i += 1 
    
  path_to_integrated_output = path_to_integ + wildcards[0] + '/integrated_cas.rds'
  path_to_combined_sce = path_to_integ + wildcards[0] + '/uncorrected_cas.rds'
  
  return sce_cas_path_list, path_to_integrated_output, path_to_combined_sce

# generate report for the integrated sample

def cell_assign_multiple_report_files(wildcards, sce_cas, uncorrected, integrated_cas_path):

    sce_cas = remove_last_two_dirs(sce_cas)
    uncorrected_path = remove_last_two_dirs(uncorrected)
    integrated_cas_path = remove_last_two_dirs(integrated_cas_path)

    separate_samples = wildcards[0].split('-')

    i = 0
    sce_cas_file_paths = []

    while i < len(separate_samples):

        sce_cas_path = sce_cas + separate_samples[i] + '/sce_norm_cas.rds'
        sce_cas_file_paths.append(sce_cas_path)

        i += 1

    path_to_uncorrected = uncorrected_path + wildcards[0] + '/uncorrected.rds'
    path_to_integrated = integrated_cas_path + wildcards[0] + '/integrated_cas.rds'

    return sce_cas_file_paths, path_to_uncorrected, path_to_integrated


# def cell_assign_multiple_report_files(wildcards): 
#     
#     separate_samples = wildcards[0].split('-')
#         
#     i = 0 
#     sce_cas_file_paths = []
#     
#     while i < len(separate_samples): 
#         
#         sce_cas_path = '../data/cellassign/' + separate_samples[i] + '/sce_norm_cas.rds'
#         sce_cas_file_paths.append(sce_cas_path)
#         
#         i += 1 
#         
#     path_to_uncorrected = '../data/integrated/' + wildcards[0] + '/uncorrected.rds'
#     path_to_integrated = '../data/integrated_cellassign/' + wildcards[0] + '/integrated_cas.rds'
#     
#     return sce_cas_file_paths, path_to_uncorrected, path_to_integrated


  
  
  













