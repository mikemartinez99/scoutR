#' @title trackCells
#' 
#' @description Generate a data-frame showing how many cells you lose for each specific parameter threshold. This function takes a named list of `Seurat` objects and a named list of expressions to pass to `seurat::subset`. This function DOES NOT directly apply the filtering. It is meant to be used a diagnostic. 
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import Signac
#' @import rlang
#' 
#' @export
#'
#' @param sampleList A named list of `Seurat` objects with valid meta.data slots
#' @param thresholds A named list of expressions to pass to `seurat::subset` (the names will serve as column names in the output dataframe)
#'
#' @returns A dataframe
#' 
#' @examples # Generate individual plots as well as a combined one
#' @examples atacList <- list("obj1", "obj2", "obj3")
#' @examples thresholds <- list(
#' @examples              pct_reads_in_peaks = "pct_reads_in_peaks > 50",
#' @examples              blacklist_fraction = "blacklist_fraction < 0.03",
#' @examples              TSS_enrichment_low = "TSS.enrichment > 2",
#' @examples              TSS_enrichment_high = "TSS.enrichment < 15",
#' @examples              nucleosome_signal = "nucleosome_signal < 4")
#' @examples # Run the function
#' @examples cellCounts <- trackCells(atacList, thresholds)

trackCells <- function(sampleList, thresholds) {
  if (!is.list(sampleList) || is.null(names(sampleList))) {
    stop("`sampleList` must be a *named* list of objects with @meta.data.")
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
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
      next
    }
    
    nCellsPre <- nrow(x@meta.data)
    removal_counts <- numeric(length(thresholds))
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
      
      removal_counts[filter_name] <- nCellsPre - nrow(x_filt@meta.data)
    }
    
    total_lost <- sum(removal_counts, na.rm = TRUE)
    
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
