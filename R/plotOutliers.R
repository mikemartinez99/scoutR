#' @title plotOutliers
#' 
#' @description Generates publication quality violin plots showing outlier 
#' distributions for any QC metric. Can be used in conjunction
#' with `ggpubr::ggarrange` to plot individual outlier plots as well as a paneled plot. 
#'
#' @param metadata A dataframe containing all metadata of interest for samples
#' @param variable A valid column in the metadata dataframe
#' @param outlierData A column name indicating results of `scuttle::isOutlier` function
#' @param outDir A directory path ending in "/" to hold output files
#' @param width Plot width
#' @param height Plot height
#' 
#' @returns A ggplot object
#' 
#' @examples # Run after using `markOutliers`
#' @examples plotOutliers(metadata = metadata, 
#' variable = "nCount_peaks", 
#' outlierData = "nCount_peaks_outliers", 
#' outDir = outputDir, 
#' width = 10, 
#' height = 10)
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import Signac
#' @import ggpubr
#' @import ggplot2
#' @import ggExtra
#' @import rlang
#' @import scuttle
#' 
#' @export


plotOutliers <- function(metadata, variable, outlierData, outDir, width, height) {
  # Convert `variable` to a symbol
  var <- ensym(variable)
  outlierVar <- ensym(outlierData)
  label <- gsub("_", " ", variable)
  
  p1 <- ggplot2::ggplot(metadata, aes(x = orig.ident, y = .data[[as_string(var)]])) +
    geom_jitter(aes(color = .data[[as_string(outlierVar)]]), width = 0.2, alpha = 0.5, size = 1.5) +
    geom_violin(fill = "gray90", color = "gray40") +
    geom_boxplot(width = 0.2, outliers = FALSE) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    theme_classic(base_size = 16) +
    labs(x = "Sample", 
         y = label, 
         color = "Outlier",
         title = label) +
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(angle = 90, face = "bold"))
  
  #----- Filename
  fileName <- paste0(var, "_outliers.png")
  ggplot2::ggsave(paste0(outDir, fileName), p1, width = width, height = height)
  
  return(p1)
  
}