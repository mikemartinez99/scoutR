#' @name trackCells
#' 
#' @title trackCells
#' 
#' @description Generate a data-frame showing how many cells you lose for each 
#' specific parameter threshold. This function takes a named list of `Seurat` 
#' objects and a named list of expressions to pass to `seurat::subset`. 
#' This function DOES NOT directly apply the filtering, It is meant to be used a diagnostic.
#' This function takes into account barcodes that are "double-dippers" across different
#' metric (i.e., a barcode that is filtered independently across multiple metrics.) 
#' 
#' @param sampleList A named list of `Seurat` objects with valid meta.data slots
#' @param thresholds A named list of expressions to pass to `seurat::subset` (the names will serve as column names in the output dataframe)
#'
#' @returns A dataframe
#' 
#' @examples # Generate individual plots as well as a combined one
#' atacList <- list("obj1", "obj2", "obj3")
#' thresholds <- list(
#'   pct_reads_in_peaks = "pct_reads_in_peaks > 50",
#'   blacklist_fraction = "blacklist_fraction < 0.03",
#'   TSS_enrichment_low = "TSS.enrichment > 2",
#'   TSS_enrichment_high = "TSS.enrichment < 15",
#'   nucleosome_signal = "nucleosome_signal < 4")
#' # Run the function
#' cellCounts <- trackCells(
#'   atacList, 
#'   thresholds)
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import Signac
#' @import rlang
#' 
#' @export

trackCells <- function(sampleList, thresholds) {
  if (!is.list(sampleList) || is.null(names(sampleList))) {
    stop("`sampleList` must be a *named* list of Seurat objects.")
  }
  
  if (!is.list(thresholds) || length(thresholds) == 0) {
    stop("`thresholds` must be a non-empty named list of filter expressions.")
  }
  
  result_cols <- c("Sample", "Original_Count", names(thresholds), "Total_Cells_Removed")
  atacFiltData <- list()
  
  for (i in names(sampleList)) {
    message("#-------------------------------------#")
    message("Tracking cells for sample ", i)
    
    x <- sampleList[[i]]
    
    if (!"meta.data" %in% slotNames(x)) {
      warning("Sample ", i, " does not have a @meta.data slot. Skipping.")
      next
    }
    
    nCellsPre <- nrow(x@meta.data)
    removal_counts <- numeric(length(thresholds))
    failed_barcodes_all <- c()
    names(removal_counts) <- names(thresholds)
    
    for (filter_name in names(thresholds)) {
      filter_expr <- thresholds[[filter_name]]
      
      parsed_expr <- tryCatch({
        rlang::parse_expr(filter_expr)
      }, error = function(e) {
        warning("Invalid filter expression for '", filter_name, "': ", filter_expr)
        return(NULL)
      })
      
      if (is.null(parsed_expr)) {
        removal_counts[filter_name] <- NA
        next
      }
      
      x_filt <- tryCatch({
        subset(x, subset = !!parsed_expr)
      }, error = function(e) {
        warning("Failed to apply filter '", filter_name, "' on sample ", i, ": ", e$message)
        return(x)
      })
      
      passed_barcodes <- colnames(x_filt)
      failed_barcodes <- setdiff(colnames(x), passed_barcodes)
      
      # Save failed barcodes
      failed_barcodes_all <- union(failed_barcodes_all, failed_barcodes)
      
      removal_counts[filter_name] <- length(failed_barcodes)
    }
    
    total_lost <- length(failed_barcodes_all)
    
    row_df <- data.frame(Sample = i,
                         Original_Count = nCellsPre,
                         as.list(removal_counts),
                         Total_Cells_Removed = total_lost,
                         stringsAsFactors = FALSE)
    
    atacFiltData[[i]] <- row_df
  }
  
  atacFiltData <- do.call(rbind, atacFiltData)
  numeric_cols <- setdiff(colnames(atacFiltData), "Sample")
  atacFiltData[numeric_cols] <- lapply(atacFiltData[numeric_cols], as.numeric)
  
  return(atacFiltData)
}
