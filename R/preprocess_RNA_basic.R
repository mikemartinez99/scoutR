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
preprocess_RNA_basic <- function(x, ident, nCount_thresh = 300, mito_gene_prefix) {
  message(ident)
  message("#------------------------------------------------------------------#")
  
  #----- Check for presence of nCount_RNA
  if ("nCount_RNA" %in% colnames(x@meta.data)) {
    
    # Count how many cells are below the threshold
    numCellsFail <- length(x@meta.data$nCount_RNA < nCount_thresh)
    
    if (numCellsFail == 0) {
      stop("No cells found")
    } else {
      message(paste0(numCellsFail, " cells were below nCount_thresh >= ", nCount_thresh))
    }
    
    # Subset based on nCount_RNA
    xSub <- subset(x, subset = nCount_RNA >= nCount_thresh)
    
  } else {
    stop("None of the requested variables were found")
  }
  
  #----- Check for mitochondrial genes
  mito_genes <- grep(mito_gene_prefix, rownames(xSub), value = TRUE)
  if (length(mito_genes) == 0) {
    warning(paste0("No mitochondrial genes found using pattern '", mito_gene_prefix,
                   "' in sample ", ident))
  } else {
    message(paste0(length(mito_genes), " mitochondrial genes found in sample ", ident))
  }
  
  #----- Add percent mitochondria column
  xSub$percentMT <- PercentageFeatureSet(object = xSub, pattern = mito_gene_prefix)
  
  #----- Add log10 genes per UMI
  xSub$log10GenesPerUMI <- log10(xSub$nFeature_RNA) / log10(xSub$nCount_RNA)
  
  #----- End sample
  message("Done")
  message("\n")
  return(xSub)
}

