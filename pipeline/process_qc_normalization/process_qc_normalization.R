# source some packages and general functions
print("Sourcing files.")

library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

# source the files that have the functions for make_sce, qc and normalization steps 
source(here('pipeline', 'process_qc_normalization', 'data_preparation', 'convert_to_sces.R')) # make sce
source(here('pipeline', 'process_qc_normalization', 'qc', 'default_qc.R')) # QC
source(here('pipeline', 'process_qc_normalization', 'normalization', 'normalize_sce.R')) # normalize

parser <- ArgumentParser(description = "Make sce, QC and normalize")

# convert to sce
parser$add_argument('--path_to_10X', metavar='DIRECTORY', type='character', help="Path to cellranger output directory")
parser$add_argument('--id', metavar='VARIABLE', type='character', help="unique id of the input dataset")

# QC
parser$add_argument('--whichMethod', metavar = 'FILE', type = 'character', help = 'pick method for qc: default or quantile')
parser$add_argument('--mito_thresh_max', metavar='FILE', type='integer', help="Maximum pct of mito counts allowed in cells")
parser$add_argument('--mito_thresh_min', metavar='FILE', type='integer', help="Minium pct of mito counts allowed in cells")
parser$add_argument('--ribo_thresh_max', metavar='FILE', type='integer', help="Maximum pct of ribo counts allowed in cells")
parser$add_argument('--nmads', metavar='FILE', type='integer', help="MAD threshold")
parser$add_argument('--min_features', metavar='FILE', type='integer', help="Minimum number of detected genes")
parser$add_argument('--remove_mito_and_ribo', metavar = 'FILE', type = 'character', help = 'Do you want to remove mito and ribo genes from the sces?')

# normalization 
parser$add_argument('--output_file_name', metavar='FILE', type='character', help="Path to normalized sce")

args <- parser$parse_args()

# convert 10X output to sce 
sce <- convert_to_sces(path_to_10X = args$path_to_10X, 
                       id = args$id)

# QC
sce <- make_sce_qc(whichMethod = args$whichMethod, 
                   sce = sce, 
                   mito_thresh_max = args$mito_thresh_max, 
                   mito_thresh_min = args$mito_thresh_min, 
                   ribo_thresh_max = args$ribo_thresh_max, 
                   nmads = args$nmads, 
                   min_features = args$min_features, 
                   remove_mito_and_ribo = args$remove_mito_and_ribo)

# normalization
sce <- normalize_sce(sce = sce, 
                     output_file_name = args$output_file_name)















