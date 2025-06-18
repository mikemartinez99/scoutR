#' @title qcViolin
#' 
#' @description Generates publication quality quality control violins. Requires that you concatenate all
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
#' @param metadata A dataframe containing all metadata of interest for samples
#' @param variable A valid column name present in your metadata dataframe
#' @param logTransform Logical, whether or not to log transform the y-axis
#' @param sampleColors A color vector for samples
#' @param figDir A directory path ending in "/" to hold output files
#'
#' @returns A ggplot object
#' 
#' @examples # Generate individual plots as well as a combined one
#' @examples arrangedViolins <- ggarrange(customViolin2(metadata, "pct_reads_in_peaks", logTransform = FALSE, colors, figDir),
#' @examples qcViolin(metadata, "atac_peak_region_fragments", logTransform = TRUE, colors, figDir),
#' @examples qcViolin(metadata, "TSS.enrichment", logTransform = FALSE, colors, figDir),
#' @examples qcViolin(metadata, "nucleosome_signal", logTransform = FALSE, colors, figDir),
#' @examples qcViolin(metadata, "nFeature_peaks", logTransform = FALSE, colors, figDir))
#' @examples ggsave(paste0(figDir, "ATAC_Violin_Summary.png"), arrangedViolins, width = 16, height = 10)


qcViolin <- function(metadata, variable, logTransform, sampleColors, figDir) {
  
  # Convert `variable` to a symbol
  variable <- ensym(variable)
  label <- gsub("_", " ", variable)
  
  
  if (logTransform == TRUE) {
    p1 <- ggplot(combinedMeta, aes(x = orig.ident, y = log(.data[[as_string(variable)]]), fill = orig.ident)) +
      geom_point(position = position_jitter(height = 0.2, width = 0.2), shape = 21, size = 1, alpha = 0.6) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(outliers = FALSE, width = 0.2) +
      theme_classic(base_size = 16) +
      labs(y = paste0("Log(", variable , ")"),
           x = "", 
           fill = "Sample",
           title = label) +
      scale_fill_manual(values = colors) +  
      theme(axis.text = element_text(angle = 90, face = "bold"),
            legend.position = "none")
    
    # Set filename
    fileName <- paste0("Log_", variable, "_QC_Violins.png")
    ggsave(paste0(violinDir, fileName), p1, width = 8, height = 8)
    
  } else {
    p1 <- ggplot(combinedMeta, aes(x = orig.ident, y = .data[[as_string(variable)]], fill = orig.ident)) +
      geom_point(position = position_jitter(height = 0.2, width = 0.2), shape = 21, size = 1, alpha = 0.6) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(outliers = FALSE, width = 0.2) +
      theme_classic(base_size = 16) +
      labs(y = variable, 
           x = "", 
           fill = "Sample",
           title = label) +
      scale_fill_manual(values = colors) +  
      theme(axis.text = element_text(angle = 90, face = "bold"),
            legend.position = "none")
    
    # Set filename
    fileName <- paste0(variable, "_QC_Violins.png")
    ggsave(paste0(figDir, fileName), p1, width = 8, height = 8)
    
  }
  return(p1)
} 
