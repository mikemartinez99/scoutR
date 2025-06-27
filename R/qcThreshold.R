#' @title qcThreshold
#' 
#' @description Generates publication quality quality control violins showing thresholds. Requires that you concatenate all
#' sample metadata across multiple objects into one dataframe. Can be used in conjunction with `ggpubr::ggarrange`
#' to plot panelled figure of all QC metrics. 
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import Signac
#' @import ggpubr
#' @import ggplot2
#' @import rlang
#' 
#' @export
#'
#' @param metadata A dataframe containing all metadata of interest for samples
#' @param variable A valid column name present in your metadata dataframe
#' @param threshold A threshold for your metadata category of interest
#' @param logTransform Logical, whether or not to log transform the y-axis
#' @param sampleColors A color vector for samples
#' @param figDir A directory path ending in "/" to hold output files
#' @param width Figure width
#' @param height Figure height
#'
#' @returns A ggplot object
#' 
#' @examples # Generate individual plots as well as a combined one
#' @examples arrangedThresholds <- ggarrange(qcThreshold(metadata, "pct_reads_in_peaks", threshold = 50, logTransform = FALSE, colors, figDir, width = 8, height = 8),
#' @examples qcThreshold(metadata, "atac_peak_region_fragments", threshold = 20000, logTransform = TRUE, colors, figDir, width = 8, height = 8),
#' @examples qcThreshold(metadata, "TSS.enrichment", threshold = 2, logTransform = FALSE, colors, figDir, width = 8, height = 8),
#' @examples qcThreshold(metadata, "nucleosome_signal", threshold = 4,  logTransform = FALSE, colors, figDir, width = 8, height = 8))
#' @examples ggsave(paste0(figDir, "ATAC_Threshold_Summary.png"), arrangedThresholds, width = 16, height = 10)


qcThreshold <- function(metadata, variable, threshold, logTransform, sampleColors, figDir, width = width, height = height) {
  
  # Convert `variable` to a symbol
  variable <- ensym(variable)
  label <- gsub("_", " ", variable)
  
  
  if (logTransform == TRUE) {
    p1 <- ggplot(metadata, aes(x = orig.ident, y = log(.data[[as_string(variable)]]), fill = orig.ident)) +
      geom_point(position = position_jitter(height = 0.2, width = 0.2), shape = 21, size = 1, alpha = 0.6) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(outliers = FALSE, width = 0.2) +
      theme_classic(base_size = 16) +
      labs(y = paste0("Log(", variable , ")"),
           x = "", 
           fill = "Sample",
           title = label) +
      geom_hline(yintercept = threshold, linetype = "dashed") +
      scale_fill_manual(values = sampleColors) +  
      theme(axis.text = element_text(angle = 90, face = "bold"),
            legend.position = "none")
    
    # Set filename
    fileName <- paste0("Log_", variable, "_QC_Threshold.png")
    ggsave(paste0(figDir, fileName), p1, width = 8, height = 8)
    
  } else {
    p1 <- ggplot(metadata, aes(x = orig.ident, y = .data[[as_string(variable)]], fill = orig.ident)) +
      geom_point(position = position_jitter(height = 0.2, width = 0.2), shape = 21, size = 1, alpha = 0.6) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(outliers = FALSE, width = 0.2) +
      theme_classic(base_size = 16) +
      labs(y = variable, 
           x = "", 
           fill = "Sample",
           title = label) +
      scale_fill_manual(values = sampleColors) +  
      geom_hline(yintercept = threshold, linetype = "dashed") +
      theme(axis.text = element_text(angle = 90, face = "bold"),
            legend.position = "none")
    
    # Set filename
    fileName <- paste0(variable, "_QC_Threshold.png")
    ggsave(paste0(figDir, fileName), p1, width = width, height = height)
    
  }
  return(p1)
} 
