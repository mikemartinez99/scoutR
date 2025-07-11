#' @name calcMAD
#' Calculate MAD-based Thresholds for a Metadata Variable
#' 
#' @title calcMAD
#' 
#' @description
#' Calculates median average deviation for a set value for a variable of 
#' interest to help quantitatively choose thresholds
#'
#' @param x A Seurat object
#' @param variable A valid column name in metadata to be tested
#' @param nmads The number of median absolute deviations to consider
#' @param bound One of "higher", or "lower"
#'
#' @returns A single numeric value representing the value at the specified nMAD
#' 
#' @examples
#' # Calculate an upper bound on a single Seurat object
#' calcMAD(
#'   x = seurat_obj,
#'   variable = "log_atac_peak_region_fragments",
#'   nmads = 5,
#'   bound = "higher")
#' 
#' @importFrom purrr map 
#' @importFrom stats mad median
#' @import Seurat
#' @import scuttle
#' 
#' @export


calcMAD <- function(x, variable, nmads, bound = c("higher", "lower")) {
  #----- Match arg for safety
  bound <- match.arg(bound)
  
  #----- Check inputs
  if (!inherits(x, "Seurat")){
    stop("`x` must be a Seurat object.")
  } 
  
  if (missing(variable)) {
    stop("You must specify a `variable` to compute MAD.")
  }
  
  if (!is.numeric(nmads) || length(nmads) != 1 || nmads < 0) {
    stop("`nmads` must be a single non-negative numeric value.")
  }
  
  #----- Get variable safely
  variable <- rlang::ensym(variable)
  var_name <- rlang::as_string(variable)
  
  meta <- x@meta.data
  
  if (!(var_name %in% colnames(meta))) {
    stop(paste0("`", var_name, "` not found in `@meta.data`."))
  }
  
  values <- meta[[var_name]]
  if (!is.numeric(values)) {
    stop(paste0("Variable `", var_name, "` must be numeric."))
  }
  
  #----- Get median and MAD
  med <- median(values, na.rm = TRUE)
  mad_val <- mad(values, constant = 1, na.rm = TRUE)
  
  if (is.na(med) || is.na(mad_val)) {
    stop("Median or MAD could not be calculated (possible NA values in data).")
  }
  
  #----- Compute bound
  if (bound == "higher") {
    result <- med - nmads * mad_val
    message("Upper Bound: ", result)
  } else {
    result <- med + nmads * mad_val
    message("Lower Bound: ", result)
  }
  
  return(result)
}