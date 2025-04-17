#' preprocess_RNA_basic
#' 
#' @description
#'Filters a seurat object based on nCount RNA and adds additional metadata for percent mitochondria and log10 genes per UMIs
#' 
#'
#' @import Seurat
#' 
#' @param x A seurat object
#' @param ident Ident corresponding to seurat object
#' @param nCount_thresh Threshold for filtering based on nCount_RNA (default = 300)
#' @param mito_gene_prefix `pattern` argument from `Seurat::PercentageFeatureSet` (i.e., "^mt-")
#'
#' @returns xSub, a seurat object filtered for `nCount_thresh` cells with added metadata for percent mitochondria and log10 genes per UMI
#' @export
#'
#' @examples # Run on a single seurat object
#' @examples subset <- preprocess_RNA_basic(x, "Sample_1", nCount_thresh = 250, mito_gene_prefix = "^mt-")
preprocess_RNA_basic <- function(x, ident, nCount_thresh=300, mito_gene_prefix) {
  print(ident)
  
  # Subset based on nCount RNA
  xSub <- subset(x,
                 subset = nCount_RNA >= nCount_thresh)
  
  
  # Add mito qc column ("^mt-")
  xSub$percentMT <- PercentageFeatureSet(object = xSub,
                                         pattern = mito_gene_prefix)
  
  # Add log10 gene/UMI (complexity measure of dataset)
  xSub$log10GenesPerUMI <- log10(xSub$nFeature_RNA) / log10(xSub$nCount_RNA)
  
  return(xSub)
}

