#' @title qcRidges
#' 
#' @description Generates publication quality quality control ridgeplots. Requires that you concatenate all
#' sample metadata across multiple objects into one dataframe. Can be used in conjunction with `ggpubr::ggarrange`
#' to plot panelled figure of all QC metrics. 
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import Signac
#' @import ggpubr
#' @import ggplot2
#' @import rlang
#' @import ggridges
#'
#' @param metadata A dataframe containing all metadata of interest for samples
#' @param variable A valid column name present in your metadata dataframe
#' @param sampleColors A color vector for samples
#' @param figDir A directory path ending in "/" to hold output files
#'
#' @returns A ggplot object
#' 
#' @examples # Generate individual plots as well as a combined one
#' @examples arrangedViolins <- ggarrange(qcRidges(metadata, "pct_reads_in_peaks", logTransform = FALSE, colors, figDir),
#' @examples qcRidges(metadata, "atac_peak_region_fragments", colors, figDir, width = 8, height = 10),
#' @examples qcRidges(metadata, "TSS.enrichment", colors, figDir, width = 8, height = 10),
#' @examples qcRidges(metadata, "nucleosome_signal", colors, figDir, width = 8, height = 10),
#' @examples qcRidges(metadata, "nFeature_peaks", colors, figDir, width = 8, height = 10)
#' @examples ggsave(paste0(figDir, "ATAC_Violin_Summary.png"), arrangedViolins, width = 16, height = 10)

qcRidges <- function(metadata, variable, logTransform, sampleColors, figDir, width, height) {
  
  # Convert `variable` to a symbol
  variable <- ensym(variable)
  label <- gsub("_", " ", variable)
  
  p1 <- ggplot(metadata, aes(x = .data[[as_string(variable)]], y = orig.ident, fill = orig.ident)) +
    geom_density_ridges(scale = 1, alpha = 0.8) +  
    theme_classic(base_size = 16) +  
    labs(x = "label", 
         y = "") +
    scale_fill_manual(values = colors) +  
    theme(axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, face = "bold"),
          legend.position = "none")
  
  # Set filename
  fileName <- paste0("Log_", variable, "_QC_Violins.png")
  ggsave(paste0(figDir, fileName), p1, width = width, height = height)
  
 
  return(p1)
} 
