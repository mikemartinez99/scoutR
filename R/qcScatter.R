#' @title qcScatter
#' 
#' @description Generates publication quality quality control scatter plot for 2 metadata variables. Requires that you concatenate all
#' sample metadata across multiple objects into one dataframe. 
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import Signac
#' @import ggpubr
#' @import ggplot2
#' @import ggExtra
#' @import rlang
#' @export
#'
#' @param metadata A dataframe containing all metadata of interest for samples
#' @param XVar A valid column name present in your metadata dataframe to plot on the X-axis
#' @param YVar A valid column name present in your metadata dataframe to plot on the Y-axis
#' @param logTransformX Logical, whether or not to log transform the X-axis (default = FALSE)
#' @param logTransformY Logical, whether or not to log transform the Y-axis (default = FALSE)
#' @param colors A color vector for samples
#' @param outDir A directory path ending in "/" to hold output files
#' @param width Plot width
#' @param height Plot height
#'
#' @returns A ggplot object
#' 
#' @examples # Generate scatter plot
#' @examples compareDists(metadata, Xvar = "log_atac_peak_region_fragments", Yvar = "pct_reads_in_peaks", logTransformX = FALSE, 
#' @examples logTransformY = FALSE, colors = colors, outDir = figDir, width = 12, height = 10)

qcScatter <- function(metadata, Xvar, Yvar, 
                         logTransformX = FALSE, logTransformY = FALSE, 
                         colors, outDir, width, height) {
  
  # Convert `variables` to a symbol
  variableX <- ensym(Xvar)
  variableY <- ensym(Yvar)
  
  # Handle dynamic axis labels
  x_label <- if (logTransformX) paste0("Log(", gsub("_", " ", as_string(variableX)), ")") else gsub("_", " ", as_string(variableX))
  y_label <- if (logTransformY) paste0("Log(", gsub("_", " ", as_string(variableY)), ")") else gsub("_", " ", as_string(variableY))
  title <- paste0(x_label, " vs ", y_label)
  
  # If Y var is logged
  p1 <- ggplot(metadata, aes(x = .data[[as_string(variableX)]], y = .data[[as_string(variableY)]], color = orig.ident)) +
    geom_point(size = 2, alpha = 0.3) +  
    theme_classic(base_size = 16) +  
    labs(x = x_label, 
         y = y_label, 
         color = "Sample",
         title = title) +
    scale_color_manual(values = colors) +  
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"),
          legend.position = "none")
  
  # Conditionally apply log scales
  if (logTransformX) p1 <- p1 + scale_x_log10()
  if (logTransformY) p1 <- p1 + scale_y_log10()
  
  #----- Add margin annotations of distributions
  p2 <- ggMarginal(p1, groupColour = TRUE, groupFill = FALSE)
  
  p3 <- ggplot(metadata, aes(x = .data[[as_string(variableX)]], y = .data[[as_string(variableY)]], color = orig.ident)) +
    geom_point(size = 2, alpha = 0.3) +  
    theme_classic(base_size = 16) +  
    labs(x = x_label, 
         y = y_label, 
         color = "Sample") +
    scale_color_manual(values = colors) +  
    facet_grid(~orig.ident) +
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "none")
  
  # Apply log scales to faceted plot if needed
  if (logTransformX) p3 <- p3 + scale_x_log10()
  if (logTransformY) p3 <- p3 + scale_y_log10()
  
  final <- cowplot::plot_grid(p2, p3, ncol = 1)
  fileName <- paste0(Xvar, "_vs_", Yvar, ".png")
  ggsave(paste0(outDir, fileName), final, width = width, height = height)
  return(final)
}

