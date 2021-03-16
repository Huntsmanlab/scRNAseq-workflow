# WHAT DOES THIS SCRIPT DO? ####
# this script only has functions that do the following: 
# make sce objects
# compute QC metrics on these sce objects
# perform QC by eliminating low quality cells 
# make plots by comparing various data sets to one another, we have another script compare data sets among themselves rather than comparing them with each other 


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ BEGIN @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# PERFORM QC FUNCTION ####
# This function perform QC and drops the low quality cells. it needs a list of sce objects to as an input.
# again, number of sce objects in the list does't matter. 
# WHAT ARE THE INPUTS OF THIS FUNCTION? 
# sce_list : list of sce objects that we ant to perform QC on 
# QCmethod : here, we decide how to perform QC. we have two options for QC method, both are strings:
# 1. based on thresholds. to pick this option: QCmethod = 'threshold'
# 2. based on outliers top pick this option: QCmethod = 'outliers'
# 3. based on 25 % quantiles. to pick this option: QCmethod == 'quantile'
# the function also has the default values for QC with thresholds. If we want to use quantiles or outliers and not thresholds, we dont need to pass these threshold values.
# but if we use the threshold methods, we may need to assign new values to these thresholds given in the function. wow much such a flexible function. 
# EXPL: perform(sces, 'quantile', libSize = 2000)

performQC <- function(sce_list, 
                      QCmethod, 
                      libSize = 1000, 
                      detectedSize = 500,
                      mitoMaxThreshold = 10, 
                      mitoMinThreshold = 2, 
                      log10GenesPerUMI_minThreshold = 0.8) {
  
  # find how many items we passed into the function. 
  max <- length(sce_list)
  
  # make empty lists to save the results of our while loop 
  filteredSCElist = list()
  filteredQCmetricslist = list() 
  
  # initiate the while loop, our first index is 1. 
  index <- 1
  
  while(index <= max) { 
    
    # find which element in the list we work on based on the index, we will perform QC on that sce. 
    sce <- sce_list[[index]]
    
    # calculate the QC metrics here on the certain sce in the list. this is outside of the if loop because regardless of the method we pick, the way we calculate 
    # the QC metrics is the same. how we treat them is different though. 
    QCresults_sce <- perCellQCMetrics(sce, subsets = list( mito = grepl("^MT-", rowData(sce)$Symbol), ribo = grepl("^RP[LS]", rowData(sce)$Symbol)))  
    
    # add one more important metric to our QC results, log10(ratio of detected genes / RNA reads)
    QCresults_sce$log10GenesPerUMI <- log10(QCresults_sce$detected) / log10(QCresults_sce$sum)
    
    # we also need to generate a data frame and apply our indices to that - we will need that to make our QC plots later on. 
    # we will update this data frame as well as we do QC. 
    filteredQCresults <- as.data.frame(QCresults_sce) # turn the S$ object into a data frame 
    
    if (QCmethod == 'threshold') {  
      
      # now do the filtering and get the row numbers according to the conditions. this gives us row numbers to keep. 
      index_keep <- which(QCresults_sce$sum > libSize &
                            QCresults_sce$detected > detectedSize &
                            QCresults_sce$subsets_mito_percent < mitoMaxThreshold & 
                            QCresults_sce$subsets_mito_percent > mitoMinThreshold &
                            QCresults_sce$log10GenesPerUMI > log10GenesPerUMI_minThreshold) 
      
      # either method creates indices for us, now we apply the indices to sce object - only keep the high quality rows. 
      # UPDATED SCE: 
      filtered_sce <- sce[index_keep, ]
      
      # updated QC metrics data frame. 
      filteredQCresults <- filteredQCresults[index_keep, ] # apply the indices to the data frame as well. 
      
    } # end of if 
    
    # this next part is for using outliers to do QC, so QCmethod == 'quantile'. we will use a special function from scater to do that. 
    else if (QCmethod == 'outliers') {
      
      # we are using quickPerCellQC() function. it calls isOutlier() to the inputs we give it to. 
      
      # lets get a summary of the reasons for discarding cells 
      # "reasons" is table full of TRUE and FALSE values. TRUE means discard, FALSE means dont discard. 
      reasons <- quickPerCellQC(QCresults_sce, percent_subsets= c("subsets_mito_percent", "log10GenesPerUMI", "subsets_ribo_percent")) 
      nDiscarded <- colSums(as.matrix(reasons)) # total number of removed cells for various reasons
      
      # Keeping the columns we DON'T want to discard - apply the results of QC to the sce object itself.
      # UPDATED SCE: 
      filtered_sce <- sce[,!reasons$discard]
      
      # get the indices of the rows we keep, these are shown as FALSE on the discard column we just made in reasons table. 
      index_keep <- which(reasons$discard == FALSE)
      
      # apply the indices to the data frame as well, meaning keep the high quality rows. 
      filteredQCresults <- filteredQCresults[index_keep, ] 
      
    } # end of else if 
    
    else { # QCmethod == 'quantile'
      
      # this method will calculate the quantiles and remove the cells if they are outside of the quantile 
      libsize_quantile <- quantile(QCresults_sce$sum, probs = seq(0, 1, 0.25)) # this tests the number of RNA reads
      gene_quantile <- quantile(QCresults_sce$detected, probs = seq(0, 1, 0.25)) # this tests the number of detected genes 
      mito_quantile <- quantile(QCresults_sce$subsets_mito_percent, probs = seq(0, 1, 0.25)) # this test the mito RNA 
      
      # determine which cells are in the 25% quantile
      libsize_drop <- QCresults_sce$sum < libsize_quantile[['25%']]
      gene_drop <- QCresults_sce$detected < gene_quantile[['25%']]
      mito_drop <- QCresults_sce$subsets_mito_percent < mito_quantile[['25%']]
      
      idx_dropped <- libsize_drop | gene_drop | mito_drop # combine the two reasons for dropping
      # idx_dropped <- which(idx_dropped == TRUE) # get the indices where we have TRUE, that corresponds to rows to drop
      
      # drop those cells in the 25% quantile (or keep the cells not in these indices ). 
      # in other words, find the indices where we have FALSE and apply that to the sce. 
      filtered_sce <- sce[!idx_dropped, ]
      
      # now deal with the QC metrics data frame, we need to drop the low quality cells from the data frame as well. 
      filteredQCresults <- filteredQCresults[!idx_dropped, ]
      
    } # end of else - END OF THE WHOLE IF LOOP 
    
    # append our updated sces and QC data frames in a lists when the if loop is over.
    filteredSCElist[[index]] <- filtered_sce
    filteredQCmetricslist[[index]] <- filteredQCresults
    
    index = index + 1 # go to the next sce in our list and repeat the steps above 
    
  } # end of while loop 
  
  # put the two things we need to return in a list to return 
  results = list(filteredSCElist, filteredQCmetricslist)
  
  return(results)  # our output is the updated list of sces and the df with only high quality cells. we need the data frame to make plots later on. 
  
} # end of function 


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# RUN THE FUNCTIONS ####
# lists that contain the data sets we need are defined at the beginning of this script 

# # make the sce objects 
# sce_list <- make_new_sce(dataset_filepath_list)
# 
# # compute the QC metrics 
# QCmetrics_list <- computeQCmetrics(sce_list)
# 
# # perform QC on sce objects 
# QCdataframes <- performQC(sce_list)

# PREPARE DATA FUNCTION ####
# USE THIS TO PREPARE THE DATA BEFORE QC ONLY. WE HAVE ANOTHER FUNCTION TO PREPARE DATA AFTER QC. they essentially do the same thing. 
# this function prepares our data for the plotting function. it extracts the columns we are intersted in from the QC data frames. 
# the inputs are as follows: 
# dfList: a list containing data frames. these dataframes are the result of QC
# colname: choose one of the three: 'sum' 'detected' 'mito'
# average: TRUE or FALSE. TRUE gives you the average value per cell, FALSE gives you the total sum across all cells. 


# 
prepare_before_data <- function(dfList, colname, average) {
  
  # find how many items we passed into the function. 
  max <- length(dfList)
  
  # make an empty list to save the results of our while loop 
  results = list()
  
  # initiate the while loop, our first index is 1. 
  index <- 1
  
  while (index <= max) {
    
    
    if(average == FALSE) { 
      
      if(colname == 'sum') {results[[index]] <- as.numeric(colSums(as.data.frame(dfList[[index]]$sum)))}
      else if(colname == 'detected') {results[[index]] <- as.numeric(colSums(as.data.frame(dfList[[index]]$detected)))}
      else if(colname == 'mito') {stop('Oh no! Please pass average == TRUE for a mito graph.')}
      
    }
    
    else {
      
      if(colname == 'sum') {results[[index]] <- as.numeric(mean(dfList[[index]]$sum))}
      else if(colname == 'detected') {results[[index]] <- as.numeric(mean(dfList[[index]]$detected))}
      else if(colname == 'mito') {results[[index]] <- as.numeric(mean(dfList[[index]]$subsets_mito_percent))}
      
    }
    
    index = index + 1
    
  } # end of while loop 
  
  return(results)
  
}

prepare_after_data <- function(dfList, colname, average) {
  
  # find how many items we passed into the function. 
  max <- length(dfList)
  
  # make an empty list to save the results of our while loop 
  results = list()
  
  # initiate the while loop, our first index is 1. 
  index <- 1
  
  while (index <= max) {
    
    if(max == 1){ # use this when we pass a list of data frames 
      
      if(average == FALSE) { 
        
        if(colname == 'sum') {results[[index]] <- as.numeric(colSums(as.data.frame(dfList[[2]][[index]]$sum)))}
        else if(colname == 'detected') {results[[index]] <- as.numeric(colSums(as.data.frame(dfList[[2]][[index]]$detected)))}
        else if(colname == 'mito') {stop('Oh no! Please pass average == TRUE for a mito graph.')}
        
      }
      
      else {
        
        if(colname == 'sum') {results[[index]] <- as.numeric(mean(dfList[[2]][[index]]$sum))}
        else if(colname == 'detected') {results[[index]] <- as.numeric(mean(dfList[[2]][[index]]$detected))}
        else if(colname == 'mito') {results[[index]] <- as.numeric(mean(dfList[[2]][[index]]$subsets_mito_percent))}
        
      } # end of else 
      
    } # end of if(typeof(dflist) = 'list')
    
    else { # use this when we pass only 1 data frame as an S4 object and not a list of data frames 
      
      if(average == FALSE) { 
        
        if(colname == 'sum') {results[[index]] <- as.numeric(colSums(as.data.frame(dfList$sum)))}
        else if(colname == 'detected') {results[[index]] <- as.numeric(colSums(as.data.frame(dfList$detected)))}
        else if(colname == 'mito') {stop('Oh no! Please pass average == TRUE for a mito graph.')}
        
      }
      
      else { 
        
        if(colname == 'sum') {results[[index]] <- as.numeric(mean(dfList$sum))}
        else if(colname == 'detected') {results[[index]] <- as.numeric(mean(dfList$detected))}
        else if(colname == 'mito') {results[[index]] <- as.numeric(mean(dfList$subsets_mito_percent))}
        
      } # end of else 
      
      
    }
    
    index = index + 1
    
  } # end of while loop 
  
  return(results)
  
}


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# MAIN PLOTTING FUNCTION ####
# this function takes the sum data from multiple data sets ans makes plots to compare them. the plots we can make are: 
# RNA reads detected, genes detected, mitochondial RNA. 
# How the heck do i use this function? very easy, follow these instructions as we go through all of the inputs one by one. 
# data_list: this is the output of the preparedata function that we have above. for more details, scroll up and refer to that function. 
# name_list: vey simple, this is a list that contains the name of the datasets we want to work on. it muct be given as a list of strings 
# colname: this is what we want to make a graph of. the options are as follows, they all must be strings: 
# sum, detected, mito. WHAT DO THEY MEAN? 
# sum: all RNA detected in the experiment 
# detected: detected genes in the dataset 
# mito: mitochondrial RNA in this dataset. NEXT INPUT:
# average: this is either TRUE or FALSE. When average == FALSE< we make graphs of all counts in all cells, so we dont average per cell. 
# when average == TRUE, we compare the values per cell between graphs. 


make_plots <- function(data_list, name_list, colname, average) {
  
  # make a new data frame with names of data sets and the values to make the plots easily later on. 
  
  if(colname == 'sum') {
    
    if(average == FALSE) {
      
      # here we make a simple data frame with the values we are interested in. 
      sums <- data.frame(
        dataset = unlist(name_list), 
        values = unlist(lapply(data_list, function(x) x/1e3))) # divide each value by 1000, for easier visualization
      
      sums$values <- round(sums$values, digits = 2) # round to 2 digits 
      
      plot <- ggplot(sums, aes(x = dataset, y = values)) + 
        geom_bar(stat = 'identity', fill = "chartreuse3", color = "black") + 
        ylab("RNA reads (thousands)") + 
        ylim(0, 180) + 
        xlab("dataset") + 
        theme_amunzur
      
    }
    
    else if(average == TRUE) {
      
      sums <- data.frame(
        dataset = unlist(name_list), 
        values = unlist(lapply(data_list, function(x) x/1e3)))  # divide each value by 1000, for easier visualization
      
      sums$values <- round(sums$values, digits = 2) # round to 2 digits 
      
      plot <- ggplot(sums, aes(x = dataset, y = values)) + 
        geom_bar(stat = 'identity', fill = "chartreuse3", color = "black") + 
        ylab("RNA reads (thousands)") + 
        ylim(0, 50) + 
        xlab("dataset") + 
        theme_amunzur
      
    }
  } # end of sum plots 
  
  # detected plots begin here 
  if (colname == 'detected') {
    
    if(average == FALSE) {
      
      detected <- data.frame(
        dataset = unlist(name_list), 
        values = unlist(lapply(data_list, function(x) x/1e6))) # divide each value by 1000, for easier visualization
      
      detected$values <- round(detected$values, digits = 2) # round to 2 digits 
      
      plot <- ggplot(detected, aes(x = dataset, y = values)) + 
        geom_bar(stat = 'identity', fill = "brown3", color = "black") +
        ylab("detected genes (millions)") + 
        xlab("dataset") + 
        ylim(0, 20) + 
        theme_amunzur
      
    }
    
    else if(average == TRUE) {
      
      detected <- data.frame(
        dataset = unlist(name_list), 
        values = unlist(lapply(data_list, function(x) x/1e3))) # divide each value by 1000, for easier visualization
      
      detected$values <- round(detected$values, digits = 2) # round to 2 digits 
      
      plot <- ggplot(detected, aes(x = dataset, y = values)) + 
        geom_bar(stat = 'identity', fill = "brown3", color = "black") +
        ylab("average detected genes (thousands)") + 
        xlab("dataset") + 
        ylim(0, 10) + 
        theme_amunzur
    }
    
  } # end of detected plots 
  
  # beginning of mito plots 
  else if(colname == 'mito'){
    
    if (average == FALSE) { stop('You silly, average must be TRUE to make a mito percentage graph.') }
    
    mito <- data.frame(
      dataset = unlist(name_list), 
      values = unlist(data_list)) # divide each value by 1000, for easier visualization
    
    mito$values <- round(mito$values, digits = 2) # round to 2 digits 
    
    plot <- ggplot(mito, aes(x = dataset, y = values)) + 
      geom_bar(stat = 'identity', fill = "deeppink1", color = "black") +
      ylab("mito RNA percent") + 
      xlab("dataset") + 
      ylim(0, 20) + 
      theme_amunzur
  }
  
  return(plot)
  
} # end of function 

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# this function is for making various plots for one dataset. when we pass a list of datasets, it will make plots for each data set in the list. 
makePlotsPerDataset <- function(name_list) {
  
  index = 1 # start the counter 
  
  max = length(name_list)
  
  plots = list() # make an empty list to save our plots later. 
  
  while (index <= max) {
    
    # pick the data frame from QC metric lists according to the index we are at, and convert to a data frame. 
    metrics_beforeQC <- as.data.frame(QCmetrics_list[[index]]) # remember this function only computes QC metrics, doesnt drop any cells 
    metrics_afterQC <- as.data.frame(QCdataframes[[2]][[index]]) # this function removes low quality cells, so after QC
    
    a <- totalCountsHist(metrics_beforeQC, metrics_afterQC) # RNA reads in the cell 
    b <- geneNumberHist(metrics_beforeQC, metrics_afterQC) # number of detected genes 
    c <- UMIandGenes(metrics_beforeQC, metrics_afterQC) # map UMI and detected genes together 
    d <- showMitoPercentPlots(metrics_beforeQC, metrics_afterQC) # mito percentage plots 
    
    plots[[index]] <- grid.arrange(a, b, c, d)
    
    index = index + 1 
  }
  
  return(plots)
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# this is the default_qc function from max. the method is the same but we no longer use deprecated functions, and it is a little more flexible. 
# whichMethods can have 2 values, both must be strings: 'default' or 'quantile'. 
# default: what we do already at the moment, only difference is we dont use a deprecated function. 
# quantile: do qc without using the default method. instead use quantiles. 
# be careful when using this one, if many cells have high mito content, this one may not work. check summary stats after qc for good practice. 
# this function will give you two outputs: sce and qc metrics dataframe. you need to specify two separate paths to save them.
# this function will also name things generically, so we dont need to specify the file names. 

make_basic_sce_qc <- function(whichMethod, path_to_sce, mito_thresh_max, mito_thresh_min, nmads, seed, min_features){
  
  # set seed for reproducibility 
  set.seed(seed)
  
  # Read the raw sce 
  sce <- readRDS(path_to_sce)
  
  # Get mitochondrial genes for QC:
  mt_genes <- grepl("^MT-", rowData(sce)$Symbol)
  ribo_genes <- grepl("^RP[LS]", rowData(sce)$Symbol)
  feature_ctrls <- list(mito = rownames(sce)[mt_genes],
                        ribo = rownames(sce)[ribo_genes])
  
  # calculate the QC metrics. this is QC metrics before eliminating any low quality cells. 
  # regardless of the methods we choose to remove low quality cells, we need to do this step. 
  metrics <- perCellQCMetrics(sce, subsets = feature_ctrls)
  
  # save this - metrics before qc:  
  saveRDS(metrics, file = 'metrics.rds')
  
  if (whichMethod == 'default') {
    
    # THIS IS HOW WE DO QC: 
    # we judge mito percent based on thresholds
    # & log-library size is 3 MADs below the median log-library size 
    # & log-transformed # of expressed genes is 3 MADs below the median
    # we will compute indices of cells to DROP. these next steps have outputs of boolean values. TRUE means drop, FALSE means don't drop. 
    
    # this combines the cells we drop because of too many or too little mito percentage 
    mito_drop <- metrics$subsets_mito_percent > mitoMaxThreshold | metrics$subsets_mito_percent < mitoMinThreshold
    libsize_drop <- isOutlier(metrics$sum, nmads=nmads, type="lower", log=TRUE)
    gene_drop <- isOutlier(metrics$detected, nmads=nmads, type="lower", log=TRUE)
    
    # apply these indices to our sce, keeping the cells with indices FALSE
    sce_qc <- sce[ , !(mito_drop | libsize_drop | gene_drop)]
    
    metrics_qc <- metrics[!(mito_drop | libsize_drop | gene_drop), ]
    
    # save the metrics after removing low quality cells 
    saveRDS(metrics_qc, file = 'metrics_qc.rds')
    
  } # end of if loop - whichMethod == 'default'  
  
  else {
    # here we will drop cells based on quantiles
    
    # this method will calculate the quantiles and remove the cells if they are outside of the quantile, we work with the same 3 methods
    libsize_quantile <- quantile(metrics$sum, probs = seq(0, 1, 0.25)) # this tests the number of RNA reads
    gene_quantile <- quantile(metrics$detected, probs = seq(0, 1, 0.25)) # this tests the number of detected genes 
    mito_quantile <- quantile(metrics$subsets_mito_percent, probs = seq(0, 1, 0.25)) # this test the mito RNA 
    
    # determine which cells are in the 25% quantile
    libsize_drop <- metrics$sum < libsize_quantile[['25%']]
    gene_drop <- metrics$detected < gene_quantile[['25%']]
    mito_drop <- metrics$subsets_mito_percent < mito_quantile[['25%']]
    
    idx_dropped <- libsize_drop | gene_drop | mito_drop # combine the three reasons for dropping
    # idx_dropped <- which(idx_dropped == TRUE) # get the indices where we have TRUE, that corresponds to rows to drop. actually we dont need this because the next line does the same thing. 
    
    
    # drop those cells in the 25% quantile (or keep the cells not in these indices ). 
    # in other words, find the indices where we have FALSE and apply that to the sce. 
    sce_qc <- sce[, !idx_dropped, ]
    
    # we drop the same cells from our qc metrics data frame as well: 
    metrics_qc <- metrics[!idx_dropped, ]
    
    # save the metrics after removing low quality cells 
    saveRDS(metrics_qc, file = 'metrics_qc.rds')
    
  } # end of else - whichMethod == 'quantile'
  
  # Filter genes
  total_counts <- Matrix::rowSums(counts(sce_qc)) 
  zero_total_counts <- total_counts == 0
  rowData(sce_qc)[paste(sce_qc$id[1], "qc_pass", sep = "_")] <- !(zero_total_counts)
  
  # Quick cluster then Normalise, cluster min.size >= max pooling size in sizes, 
  # in this case, 150 > 101 so it should be ok
  # Each cluster should contain a sufficient number of cells for pooling, ensure by min.size = 200 
  # Careful definition of subclusters is not required for computeSumFactors
  # NOTE: output of running a function with subset.row will always be the same as 
  # the output of subsetting x beforehand and passing it into the function
  # gene_qc <- colnames(rowData(sce_qc))[grepl("_qc_pass" , colnames(rowData(sce_qc)))]
  # subset_row <- rowData(sce_qc)[, gene_qc]
  # set.seed(1234) # for IrlbaParam()
  # clusters <- quickCluster(sce_qc, min.size = 200, min.mean= 0.1, subset.row = subset_row,
  #                         use.ranks=FALSE, BSPARAM=IrlbaParam())
  # clusters <- quickCluster(sce_qc, assay.type="counts", min.size=150, min.mean= 0.1, method="igraph", use.ranks=TRUE)
  num_cells <- dim(sce)[2]
  min_size <- min(150, floor(dim(sce)[2] * 0.3))
  max_win <- min(101, min_size + 1)
  clusters <- quickCluster(sce_qc, assay.type="counts", min.size=min_size, min.mean= 0.1, method="igraph", use.ranks=FALSE, BSPARAM=IrlbaParam())
  sce_qc <- computeSumFactors(sce_qc, assay.type="counts", sizes=seq(21, max_win, 5), min.mean= 0.1, clusters=clusters)
  
  # Ensure that size factors are non-zero and non-negative before normalizing 
  if(summary(sizeFactors(sce_qc))[['Min.']] <= 0){
    stop('Negative size factors, cannot progess')
  }
  
  # Normalize using size factors (within batch)
  sce_qc <- logNormCounts(sce_qc)
  
  # Run UMAP, TSNE, PCA for visualization downstream steps: clustering, cellassign, ect.
  sce_qc <- runPCA(sce_qc, exprs_values = "logcounts", ncomponents = 200)
  sce_qc <- runTSNE(sce_qc, exprs_values = "logcounts", ntop = 500, ncomponents = 3)
  sce_qc <- runUMAP(sce_qc, exprs_values = "logcounts", ntop = 500, ncomponents = 3,
                    min_dist = 0.01, n_neighbors = 15, metric = "euclidean")
  
  # Save sce_qc as .rds file
  saveRDS(sce_qc, file = 'sce_qc.rds')
  
}


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

