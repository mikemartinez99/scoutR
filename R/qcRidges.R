#' @name qcRidges
#' 
#' @title qcRidges
#' 
#' @description Generates publication quality quality control ridge plots. 
#' Can be used in conjunction with `ggpubr::ggarrange`
#' to plot paneled figure of all QC metrics. 
#'
#' @param metadata A dataframe containing all metadata of interest for samples
#' @param variable A valid column name present in your metadata dataframe
#' @param sampleColors A color vector for samples
#' @param figDir A directory path ending in "/" to hold output files
#'
#' @returns A ggplot object
#' 
#' @examples # Generate individual plots as well as a combined one
#' arrangedViolins <- ggarrange(
#'   qcRidges(metadata, 
#'     "pct_reads_in_peaks", 
#'     logTransform = FALSE, 
#'     colors, 
#'     figDir,
#'     width = 8,
#'     height = 10),
#'   qcRidges(metadata, 
#'     "atac_peak_region_fragments", 
#'     colors, 
#'     figDir, 
#'     width = 8, 
#'     height = 10)
#'    )
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import Signac
#' @import ggpubr
#' @import ggplot2
#' @import rlang
#' @import ggridges
#' 
#' @export

qcRidges <- function(metadata, variable, sampleColors, figDir, width, height) {
  
  # Convert `variable` to a symbol
  variable <- ensym(variable)
  label <- gsub("_", " ", variable)
  
  p1 <- ggplot::ggplot(metadata, aes(x = .data[[as_string(variable)]], y = orig.ident, fill = orig.ident)) +
    geom_density_ridges(scale = 1, alpha = 0.8) +  
    theme_classic(base_size = 16) +  
    labs(x = "label", 
         y = "") +
    scale_fill_manual(values = sampleColors) +  
    theme(axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          legend.position = "none")
  
  # Set filename
  fileName <- paste0("Log_", variable, "_QC_Violins.png")
  ggplot::ggsave(paste0(figDir, fileName), p1, width = width, height = height)
  
 
  return(p1)
} 
