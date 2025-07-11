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
#' @param plottype One of "violin" or "ridgeplot"
#' @param width Figure width
#' @param height Figure height
#'
#' @returns A ggplot object
#' 
#' @examples # Generate individual plots as well as a combined one
#' @examples arrangedThresholds <- ggarrange(qcThreshold(metadata, "pct_reads_in_peaks", threshold = 50, logTransform = FALSE, colors, figDir, plottype = c("violin"), width = 8, height = 8),
#' @examples qcThreshold(metadata, "atac_peak_region_fragments", threshold = 20000, logTransform = TRUE, colors, figDir, plottype = c("violin"), width = 8, height = 8),
#' @examples qcThreshold(metadata, "TSS.enrichment", threshold = 2, logTransform = FALSE, colors, figDir, plottype = c("ridgeplot"), width = 8, height = 8),
#' @examples qcThreshold(metadata, "nucleosome_signal", threshold = 4,  logTransform = FALSE, colors, figDir, plottype = c("ridgeplot"), width = 8, height = 8))
#' @examples ggsave(paste0(figDir, "ATAC_Threshold_Summary.png"), arrangedThresholds, width = 16, height = 10)


qcThreshold <- function(metadata, variable, threshold, logTransform, sampleColors, figDir,
                        plottype = c("violin", "ridgeplot"), width = 8, height = 8) {
  
  # Match arg for safety
  plottype <- match.arg(plottype)
  
  # Convert variable to symbol and label
  variable <- rlang::ensym(variable)
  var_name <- rlang::as_string(variable)
  label <- gsub("_", " ", var_name)
  
  # Apply log transformation if requested
  if (logTransform) {
    metadata[[var_name]] <- log(metadata[[var_name]])
    ylabel <- paste0("Log(", var_name, ")")
    filename_prefix <- "Log_"
  } else {
    ylabel <- var_name
    filename_prefix <- ""
  }
  
  # Build plots
  if (plottype == "violin") {
    p <- ggplot(metadata, aes(x = orig.ident, y = .data[[var_name]], fill = orig.ident)) +
      geom_point(position = position_jitter(width = 0.2, height = 0.2), shape = 21, size = 1, alpha = 0.6) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(outlier.shape = NA, width = 0.2) +
      theme_classic(base_size = 16) +
      labs(y = ylabel, x = "", fill = "Sample", title = label) +
      geom_hline(yintercept = threshold, linetype = "dashed") +
      scale_fill_manual(values = sampleColors) +
      theme(axis.text.x = element_text(angle = 90, face = "bold"),
            legend.position = "none")
  } else if (plottype == "ridgeplot") {
    suppressWarnings(require(ggridges))  # ensure ggridges is loaded
    p <- ggplot(metadata, aes(x = .data[[var_name]], y = orig.ident, fill = orig.ident)) +
      geom_density_ridges(alpha = 0.7, scale = 1) +
      geom_vline(xintercept = threshold, linetype = "dashed") +
      theme_classic(base_size = 16) +
      labs(x = ylabel, y = "Sample", fill = "Sample", title = label) +
      scale_fill_manual(values = sampleColors) +
      theme(legend.position = "none")
  }
  
  # Save plot
  fileName <- paste0(filename_prefix, var_name, "_QC_", plottype, ".png")
  ggsave(file.path(figDir, fileName), p, width = width, height = height)
  
  return(p)
}
