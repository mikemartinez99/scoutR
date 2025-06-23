#' @title markOutliers
#' 
#' @description
#' Wrapper function for `scuttle::isOutlier`. If multiple samples need to be used, utilize `purrr::map`
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import scuttle
#' @export
#'
#' @param x A seurat object
#' @param vars A vector column names in seurat object metadata slot to be tested for outliers
#' @param nmad The number of median absolute deviations to consider a cell an outlier. 3-5 is recommended. 3 is more stringent. 
#' @param bound. One of "higher", "lower", or "both"
#'
#'
#' @returns A list of seurat objects (one per sample)
#'
#' @examples # Call the function on a single sample
#' @examples outlierObj <- markoutliers(seuratObj, vars = c("nCount_RNA", "nFeature_RNA", "percentMT", "log10GenesPerUMI"), nmad = 5, bound = "higher")
#'
#' @examples # Call the function across multiple Seurat objects
#' @examples OutlierList <- map(samples, ~ markOutliers(.x, vars = c("nCount_RNA", "nFeature_RNA", "percentMT", "log10GenesPerUMI"), nmad = 5, bound = "higher"))


markOutliers <- function(x, vars = c(...), nmad, bound) {
  #----- Input checks
  if (!inherits(x, "Seurat")) {
    stop("Input `x` must be a Seurat object.")
  }
  
  if (!is.numeric(nmad) || length(nmad) != 1 || nmad <= 0) {
    stop("`nmad` must be a single positive numeric value.")
  }
  
  valid_bounds <- c("lower", "upper", "both")
  if (!(bound %in% valid_bounds)) {
    stop("`bound` must be one of: 'lower', 'upper', 'both'.")
  }
  
  if (missing(vars) || length(vars) == 0) {
    stop("Please provide at least one variable name in `vars`.")
  }
  
  missing_vars <- vars[!vars %in% colnames(x@meta.data)]
  if (length(missing_vars) > 0) {
    stop("The following variables were not found in metadata: ",
         paste(missing_vars, collapse = ", "))
  }
  
  #----- Run
  for (i in vars) {
    message(paste0("Finding ", nmad, " outliers in ", i))
    resCol <- paste0(i, "_outliers")
    meta <- x@meta.data
    
    #----- Calculate outliers
    outliers <- isOutlier(meta[[i]], type = bound, nmads = nmad)
    meta[[resCol]] <- outliers
  }
  
  #----- Re-append metadata
  x@meta.data <- meta
  return(x)
}