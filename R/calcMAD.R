#' @title calcMAD
#' 
#' @description
#' Calculates median average deviation for a set value for a variable of interest to help quantitatively choose thresholds
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import scuttle
#' @export
#'
#' @param x A metadata table
#' @param vars A vector column names in metadata to be teted
#' @param nmad The number of median absolute deviations to consider.
#' @param bound One of "higher", "lower", or "both"
#'
#'
#' @returns Values per sample to the console
#' 
#' 