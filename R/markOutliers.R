#' @name markOutliers
#' Creates new outlier columns in Seurat metadata slot for variables of interest.
#' 
#' @title markOutliers
#' 
#' @description
#' Wrapper function for `scuttle::isOutlier`. Calculates outliers based on 
#' median-average-deviations (MAD). Can be parallelized across multiple samples
#' using `purrr::map`
#' 
#'
#' @param x A seurat object
#' @param vars A vector column names in Seurat object metadata slot to be tested
#' @param nmad The number of median absolute deviations to consider a cell an outlier. 
#'
#' @returns A list of seurat objects (one per sample)
#'
#' @examples
#' # Call the function on a single Seurat object
#' outlierObj <- markOutliers(
#'   seuratObj,
#'   vars = c("nCount_RNA", "nFeature_RNA", "percentMT", "log10GenesPerUMI"),
#'   nmad = 5,
#'   bound = "higher"
#' )
#'
#' # Apply the function across a list of Seurat objects using purrr::map
#' OutlierList <- purrr::map(
#'   samples,
#'   ~ markOutliers(
#'     .x,
#'     vars = c("nCount_RNA", "nFeature_RNA", "percentMT", "log10GenesPerUMI"),
#'     nmad = 5,
#'     bound = "higher"
#'   )
#' )
#'
#' @importFrom purrr map 
#' @import Seurat
#' @import scuttle
#' 
#' @export

markOutliers <- function(x, vars = c(...), nmad, bound) {
  #----- Input checks
  if (!inherits(x, "Seurat")) {
    stop("Input `x` must be a Seurat object.")
  }
  
  if (!is.numeric(nmad) || length(nmad) != 1 || nmad <= 0) {
    stop("`nmad` must be a single positive numeric value.")
  }
  
  valid_bounds <- c("lower", "higher", "both")
  if (!(bound %in% valid_bounds)) {
    stop("`bound` must be one of: 'lower', 'higher', 'both'.")
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
    message(paste0("Finding ", nmad, " nmad outliers in ", i))
    resCol <- paste0(i, "_outliers")
    meta <- x@meta.data
    
    #----- Calculate outliers
    outliers <- scuttle::isOutlier(meta[[i]], type = bound, nmads = nmad)
    meta[[resCol]] <- outliers
  }
  
  #----- Re-append metadata
  x@meta.data <- meta
  return(x)
}