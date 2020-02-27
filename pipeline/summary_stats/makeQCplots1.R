############################################################################################################################################
############################################################################################################################################
# NOTE TO SELF: 
# THIS SCRIPT IS TRASH. DONT USE IT. USE THE OTHER QC FUNCTIONS SCRIPT. 
############################################################################################################################################
############################################################################################################################################




# this script has a function that makes QC plots for 416b data set. The input of the function is the DATA FRAME that we get as a 
# result of calculating the QC metrics. We will call that metrics_beforeQC and metrics_afterQC. After means after we remove the low
# quality cells from the data frame. 

# libSize <- 1000
# detectedSize <- 500
# mitoMaxThreshold <- 10
# mitoMinThreshold <- 2
# log10GenesPerUMI_minThreshold <- 0.8
# riboMaxThreshold <- 10
# riboMinThreshold <- 2
# 
# # first, we need convert out QC control results to proper form: data frame 
# metrics_beforeQC <- as.data.frame(metrics_beforeQC)
# metrics_afterQC <- as.data.frame(metrics_afterQC)


# this function makes two histograms and displays them in the same window: library sizes of cells before and after QC. 
totalCountsHist <- function(metrics_beforeQC, metrics_afterQC) {
  
  # 1. HISTOGRAM - make a histogram to visualize library size of all cells - all transcripts found in cells. to make this plot, we need to run the perCellQCdataset function on sce first.  
  total_countsBEFORE <- ggplot(metrics_beforeQC, aes(x = sum / 1e3)) + 
    geom_histogram(fill = "chartreuse3", color = "black", bins = 20) + 
    ylab("Number of cells") + 
    xlab("library size (thousands)") + 
    xlim(0, 100) + 
    ylim(0, 1500) +  # lets make sure our scaleing is the same, this will make it easier to compare the two histograms 
    theme_amunzur
  
  total_countsAFTER <- ggplot(metrics_afterQC, aes(x = sum / 1e3)) + 
    geom_histogram(fill = "chartreuse3", color = "black", bins = 20) + 
    ylab("Number of cells") + 
    xlab("library size (thousands)") + 
    xlim(0, 100) + 
    ylim(0, 1500) +
    theme_amunzur 
  
  result <- list(total_countsBEFORE, total_countsAFTER)
  
  # now we visualize the two cases (beore and after removing low quality cells) side by side to see the difference 
  # return(gridExtra::grid.arrange(total_countsBEFORE, 
  #                                total_countsAFTER))
  
  return(result)
}

# run the function 
# totalCountsHist(metrics_beforeQC, metrics_afterQC)


# =============================================================================================================

 
# this function makes two histograms showing the gene number: before and after QC, we also visualize both in the same window. 
geneNumberHist <- function(metrics_beforeQC, metrics_afterQC) {
  
  # 2. HISTOGRAM - visualize the number of genes in cells 
  expressed_genesBEFORE <- ggplot(metrics_beforeQC, aes(x = detected / 1e3)) +  # number of genes is given under the column named "detected" 
    geom_histogram(fill = "brown3", color = "black", bins = 20) + 
    ylab("Number of cells") + 
    xlab("number of detected genes (thousands)") + 
    ylim(0, 1200) + 
    scale_x_continuous(limits = c(0, 12), breaks = scales::pretty_breaks(n = 10)) + 
    theme_amunzur 
  
  expressed_genesAFTER <- ggplot(metrics_afterQC, aes(x = detected / 1e3)) + 
    geom_histogram(fill = "brown3", color = "black", bins = 20) + 
    ylab("Number of cells") + 
    xlab("number of detected genes (thousands)") + 
    ylim(0, 1200) +
    scale_x_continuous(limits = c(0, 12), breaks = scales::pretty_breaks(n = 10)) + 
    theme_amunzur
  
  # return them in the same display window 
  return(gridExtra::grid.arrange(expressed_genesBEFORE, 
                                 expressed_genesAFTER))
}

# run the function 
# geneNumberHist(metrics_beforeQC, metrics_afterQC)


# =============================================================================================================


# 3. this function does the same thing as the function #1, but instead of a histogram, we use a density plot.
geneNumberDensity <- function(metrics_beforeQC, metrics_afterQC) {
  
  totalBEFORE <- ggplot(metrics_beforeQC, aes(x = sum / 1e6)) + 
    geom_density(size = 2, color = 'maroon') + 
    ylab("log10 cell density") +
    xlab('number of transcripts per cell (millions)') + 
    ggtitle('number of total transcripts per cell - BEFORE QC') + 
    scale_x_continuous(limit = c(0, 0.10), breaks = scales::pretty_breaks(n = 6)) + 
    theme_amunzur 
  
  totalAFTER <- ggplot(metrics_afterQC, aes(x = sum / 1e6)) + 
    geom_density(size = 2, color = 'maroon') + 
    ylab("log10 cell density") +
    xlab('number of transcripts per cell (millions)') + 
    ggtitle('number of total transcripts per cell - AFTER QC') + 
    scale_x_continuous(limit = c(0, 0.10), breaks = scales::pretty_breaks(n = 6)) + 
    theme_amunzur
  
  
  return(gridExtra::grid.arrange(totalBEFORE, 
                                 totalAFTER))
} 

# lets run our function 
# geneNumberDensity(metrics_beforeQC, metrics_afterQC)


# =============================================================================================================


# 4. same thing as function number 2: visualize the number of detected genes, but this time use a density function, not a histogram

# this is a density plot: 
geneDensityPlot <- function(metrics_beforeQC, metrics_afterQC) {
  
  geneDensityBEFORE <- ggplot(metrics_beforeQC, aes(x = detected)) + 
    geom_density(size = 2, color = 'springgreen3') + 
    ylab('density of genes') +
    xlab('number of detected genes per cell') + 
    scale_x_continuous(limit = c(0, 11000), breaks = scales::pretty_breaks(n = 6)) +
    ggtitle('genes per cell - BEFORE QC') + 
    theme_amunzur
  
  geneDensityAFTER <- ggplot(metrics_afterQC, aes(x = detected)) + 
    geom_density(size = 2, color = 'springgreen3') + 
    ylab('density of genes') +
    xlab('number of detected genes per cell') + 
    scale_x_continuous(limit = c(0, 11000), breaks = scales::pretty_breaks(n = 6)) +
    ggtitle('genes per cell - AFTER QC') + 
    theme_amunzur
  
  return(gridExtra::grid.arrange(geneDensityBEFORE, 
                                 geneDensityAFTER))
  
}

# geneDensityPlot(metrics_beforeQC, metrics_afterQC)



# =============================================================================================================


# 5. this is a box plot of the same parameters as the function above, showing detected genes 

genesBoxPlot <- function(metrics_beforeQC, metrics_afterQC) { 
  
  genesBoxBEFORE <- ggplot(metrics_beforeQC, aes(y = detected)) + 
    geom_boxplot() + 
    theme_amunzur + 
    ylab("number of detected genes") + 
    ggtitle("BEFORE QC")
  theme(axis.text.x = element_text(angle = 40))
  
  genesBoxAFTER <- ggplot(metrics_afterQC, aes(y = detected)) + 
    geom_boxplot() + 
    theme_amunzur + 
    ylab("number of detected genes") + 
    ggtitle("AFTER QC")
  theme(axis.text.x = element_text(angle = 40))
  
  return(gridExtra::grid.arrange(genesBoxBEFORE, 
                                 genesBoxAFTER))
  
} 

# genesBoxPlot(metrics_beforeQC, metrics_afterQC)


# =============================================================================================================


# 6. now we will visualize the ratio of detected and sum. this will help us see the ratio of genes per UMI. 
# but we will log normalize this ratio. this new column should already be added to the data frame. 

SHOWlog10GenesPerUMI <- function(metrics_beforeQC, metrics_afterQC) { 
  
  log10GenesPerUMI_BEFORE <- ggplot(metrics_beforeQC, aes(x = log10GenesPerUMI)) +
    geom_density(size = 1) + 
    ggtitle('BEFORE QC') + 
    geom_vline(xintercept = 0.80, color = 'red') + 
    scale_x_continuous(limit = c(0.7, 1.0), breaks = scales::pretty_breaks(n = 10)) +
    ylim(0, 30) + 
    theme_amunzur
  
  log10GenesPerUMI_AFTER <- ggplot(metrics_afterQC, aes(x = log10GenesPerUMI)) +
    geom_density(size = 1) + 
    ggtitle('showing log10(genes per UMI)') + 
    geom_vline(xintercept = 0.80, color = 'red') + 
    scale_x_continuous(limit = c(0.7, 1.0), breaks = scales::pretty_breaks(n = 10)) +
    ylim(0, 30) + 
    ylab('showing log10(genes per UMI)')
    theme_amunzur
  
  
  return(gridExtra::grid.arrange(log10GenesPerUMI_BEFORE, 
                                 log10GenesPerUMI_AFTER))
  
} 

# genesBoxPlot(metrics_beforeQC, metrics_afterQC)


# =============================================================================================================



# 7. compare total number of UMIs to number of genes: visualizing the correlation between 
# genes detected and number of UMIs. 
# poor cells: both low, on low end of graph
# good cells: both high, on the high end of graph 

UMIandGenes <- function(metrics) {
  
  show_both_UMI_and_genes <- ggplot(metrics, aes(y = sum / 1e3, x = detected / 1e3, color = "salmon2")) + 
    geom_point() + 
    geom_smooth(method = lm, color = "black") + 
    scale_x_continuous(limit = c(0, 15), breaks = scales::pretty_breaks(n = 10)) +
    ylab('RNA reads (thousands)') + 
    xlab('genes detected (thousands)') + 
    ylim(0, 160) + 
    theme_amunzur
  
  return(UMIandGenes)
  
}

# =============================================================================================================


# 8. mitochondrial counts: 
# this function makes two graphs (before / after QC) to show the ratio of mito RNA percent. 

showMitoPercentPlots <- function(metrics) {
  
  mito_percent_plot <- ggplot(metrics, aes(x = subsets_mito_percent)) + 
    geom_histogram(fill = 'deeppink1', color = 'black', bins = 15, position = position_nudge(x = 0.0001)) + 
    geom_vline(xintercept = 10, color = 'blue3', size = 1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + 
    xlab('mito RNA percent') + 
    ylab('number of cells') + 
    xlim(0, 35) + 
    ylim(0, 1800) + 
    ggtitle('percent of mito RNA - BEFORE QC') + 
    theme_amunzur 
  
  
  return(mito_percent_plot)
  
} 

eGFPdetectionRatePlot <- function(metrics) {
  
  mito_percent_plot <- ggplot(metrics, aes(x = subsets_mito_percent)) + 
    geom_histogram(fill = 'deeppink1', color = 'black', bins = 15, position = position_nudge(x = 0.0001)) + 
    geom_vline(xintercept = 10, color = 'blue3', size = 1) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) + 
    xlab('mito RNA percent') + 
    ylab('number of cells') + 
    xlim(0, 35) + 
    ylim(0, 1800) + 
    ggtitle('percent of mito RNA - BEFORE QC') + 
    theme_amunzur 
  
  
  return(mito_percent_plot)
  
} 


