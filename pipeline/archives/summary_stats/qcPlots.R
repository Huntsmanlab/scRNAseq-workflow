# we have some qc plots here that i used to make the summary stats. 

totalCountsPlot <- function(df, ymax) {
  
  # 1. HISTOGRAM - make a histogram to visualize library size of all cells - all transcripts found in cells. 
  total_counts <- ggplot(df, aes(x = sum / 1e3)) + 
    geom_histogram(fill = "chartreuse3", color = "black", bins = 20) + 
    geom_vline(xintercept = mean(df$sum / 1e3), color = 'red', size = 1) + 
    ylab("Number of cells") + 
    xlab("library size (thousands)") +  
    ylim(0, ymax) + 
    theme_amunzur
  
  return(total_counts) } 

geneNumberPlot <- function(df, ymax) {
  
  # 2. HISTOGRAM - visualize the number of genes in cells 
  expressed_genes <- ggplot(df, aes(x = detected / 1e3)) +  # number of genes is given under the column named "detected" 
    geom_histogram(fill = "brown3", color = "black", bins = 20) + 
    geom_vline(xintercept = mean(df$detected / 1e3), color = 'red', size = 1) + 
    ylab("Number of cells") + 
    xlab('detected genes (thousands)') + 
    ylim(0, ymax) + 
    scale_x_continuous(limits = c(0, 12), breaks = scales::pretty_breaks(n = 10)) + 
    theme_amunzur 
  
  return(expressed_genes) } 

UMIandGenesPlot <- function(df) {
  
  show_both_UMI_and_genes <- ggplot(df, aes(y = sum / 1e3, x = detected / 1e3, color = "salmon2")) + 
    geom_point() + 
    geom_smooth(method = lm, color = "black") + 
    scale_x_continuous(limit = c(0, 10), breaks = scales::pretty_breaks(n = 10)) +
    ylab('RNA reads (thousands)') + 
    xlab('genes detected (thousands)') + 
    ylim(0, 120) + 
    theme_amunzur
  
  return(show_both_UMI_and_genes) } 

MitoPercentPlot <- function(df, ymax) {
  
  mito_percent_plot <- ggplot(df, aes(x = subsets_mito_percent)) + 
    geom_histogram(fill = 'deeppink1', color = 'black', bins = 20, position = position_nudge(x = 0.0001)) + 
    geom_vline(xintercept = 10, color = 'blue3', size = 1) +
    xlab('mito RNA percent') + 
    ylab('number of cells') + 
    xlim(0, 80) + 
    ylim(0, ymax) + 
    theme_amunzur 
  
  return(mito_percent_plot) } 


################################################################################################################################################
################################################################################################################################################

