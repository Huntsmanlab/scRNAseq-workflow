# this is the default_qc function from max. the method is the same but we no longer use deprecated functions, and it is a little more flexible. 
# whichMethods can have 2 values, both must be strings: 'default' or 'quantile'. 
# default: what we do already at the moment, only difference is we dont use a deprecated function. 
# quantile: do qc without using the default method. instead use quantiles. 
# be careful when using this one, if many cells have high mito content, this one may not work. check summary stats after qc for good practice. 
# this function will give you two outputs: sce and qc metrics dataframe. you need to specify two separate paths to save them.
# this function will also name things generically, so we dont need to specify the file names. 

# source('/huntsman/amunzur/DH_organoid/pipeline/sourceFiles/utilities.R')

library(here)
source(here('pipeline', 'sourceFiles', 'utilities.R'))

parser <- ArgumentParser(description = "Do sce")
parser$add_argument('--whichMethod', metavar = 'FILE', type = 'character', help = 'pick method for qc: default or quantile')
parser$add_argument('--mito_thresh_max', metavar='FILE', type='integer', help="Maximum pct of mito counts allowed in cells")
parser$add_argument('--mito_thresh_min', metavar='FILE', type='integer', help="Minium pct of mito counts allowed in cells")
parser$add_argument('--ribo_thresh_max', metavar='FILE', type='integer', help="Maximum pct of ribo counts allowed in cells")
parser$add_argument('--nmads', metavar='FILE', type='integer', help="MAD threshold")
parser$add_argument('--min_features', metavar='FILE', type='integer', help="Minimum number of detected genes")
parser$add_argument('--remove_mito_and_ribo', metavar = 'FILE', type = 'character', help = 'Do you want to remove mito and ribo genes from the sces?')
parser$add_argument('--sce', metavar = 'DIRECTORY', type = 'character', help = 'Path to sce')
parser$add_argument('--output_file_name', metavar = 'DIRECTORY', type = 'character', help = 'Path to sce_qc')
  
args <- parser$parse_args()

print("Started QC.")
  
# set seed for reproducibility 
set.seed(1122)
  
# Read the raw sce 
sce <- readRDS(args$sce)
  
# Get mitochondrial genes for QC:
if ( grepl("^ENSMUS", sce@rowRanges@elementMetadata@listData[["ID"]][1]) ){               # For mouse
  mt_genes <- grepl("^mt", rowData(sce)$Symbol)
  ribo_genes <- grepl("^Rp[sl][[:digit:]]", rowData(sce)$Symbol)
} else {                                                                                  # For human
  mt_genes <- grepl("^MT-", rowData(sce)$Symbol)
  ribo_genes <- grepl("^RP[LS]", rowData(sce)$Symbol)
}
  
feature_ctrls <- list(mito = rownames(sce)[mt_genes],
                      ribo = rownames(sce)[ribo_genes])
  
if (args$whichMethod == 'default') {
    
    # THIS IS HOW WE DO QC: 
    # we judge mito percent based on thresholds
    # & log-library size is 3 MADs below the median log-library size 
    # & log-transformed # of expressed genes is 3 MADs below the median
    # we will compute indices of cells to DROP. these next steps have outputs of boolean values. TRUE means drop, FALSE means don't drop. 
    
    # this combines the cells we drop because of too many or too little mito percentage 
    mito_drop <- sce$subsets_mito_percent > args$mito_thresh_max | sce$subsets_mito_percent < args$mito_thresh_min
    ribo_drop <- sce$subsets_ribo_percent > args$ribo_thresh_max
    libsize_drop <- isOutlier(sce$sum, nmads=args$nmads, type="lower", log=TRUE)
    gene_drop <- isOutlier(sce$detected, nmads=args$nmads, type="lower", log=TRUE)
    
    # apply these indices to our sce, keeping the cells with indices FALSE
    sce_qc <- sce[ , !(mito_drop | libsize_drop | gene_drop | ribo_drop)]
    
  } else {
    # here we will drop cells based on quantiles
    
    # this method will calculate the quantiles and remove the cells if they are outside of the quantile, we work with the same 3 methods
    libsize_quantile <- quantile(sce$sum, probs = seq(0, 1, 0.25)) # this tests the number of RNA reads
    gene_quantile <- quantile(sce$detected, probs = seq(0, 1, 0.25)) # this tests the number of detected genes 
    mito_quantile <- quantile(sce$subsets_mito_percent, probs = seq(0, 1, 0.25)) # this test the mito RNA 
    
    # determine which cells are in the 25% quantile
    libsize_drop <- sce$sum < libsize_quantile[['25%']]
    gene_drop <- sce$detected < gene_quantile[['25%']]
    mito_drop <- sce$subsets_mito_percent < mito_quantile[['25%']]
    
    idx_dropped <- libsize_drop | gene_drop | mito_drop # combine the three reasons for dropping
    # idx_dropped <- which(idx_dropped == TRUE) # get the indices where we have TRUE, that corresponds to rows to drop. actually we dont need this because the next line does the same thing. 
    
    # drop those cells in the 25% quantile (or keep the cells not in these indices ). 
    # in other words, find the indices where we have FALSE and apply that to the sce. 
    sce_qc <- sce[, !idx_dropped, ]
    
  } 
  
  ####################################################################
  # deal with mito and ribo 
if (args$remove_mito_and_ribo == "yes") {
    
    # get index of mito genes 
    mito_idx <- which(mt_genes)
    
    # subset the sce 
    sce_qc <- sce_qc[-mito_idx, ]
    
    # get index of ribo genes 
    ribo_idx <- which(ribo_genes)
    
    # subset the sce 
    sce_qc <- sce_qc[-ribo_idx, ]  } # end of if
  
print("Finished QC.")
  
saveRDS(sce_qc, file = args$output_file_name)

