#' @title read_10x_GEX
#' 
#' @description
#' When multiple samples are present, `read_10x_GEX` can read each sample's filtered feature barcode matrix h5 file
#' and generate a seurat object for each. Each of these seurat objects is then combined into a list of objects.
#' This function can be easily parallelized using `purrr::map`
#' 
#' @importFrom purrr map 
#' @import Seurat
#' 
#' @export
#'
#' @param sampleDir Path to the sample directory containing multiple subfolders (one per sample)
#' @param sample A vector of sample names that correspond to how subfolders in `sampleDir` are named
#' @param ident A vector of desired ident names for seurat object. Idents must correspond to same order in `sample`
#' @param min.cells Minimum number of cells a feature needs to be present in to be retained
#' @param min.features Minimum number of features needed in a cell to be retained
#'
#' @returns A list of seurat objects (one per sample)
#'
#' @examples # Path to sample directory
#' @examples wd <- "/my/path/to/data/"
#' @examples # Vector of sample IDs
#' @examples sampleIDs <- c("S1", "S2", "S3")
#' @examples # Vector of idents
#' @examples idents <- c("A", "B", "C")
#' @examples # Call the function
#' @examples seuratList <- map(samples, ~ read_10x_multiGEX(.x, ident = idents, min.cells = 10, min.features = 100))

read_10x_GEX <- function(sampleDir, sample, ident, min.cells, min.features) {
  print(sample)
  counts <- Read10X_h5(filename = paste0(sampleDir, sample, "/outs/filtered_feature_bc_matrix.h5"))
  x <- CreateSeuratObject(counts = counts$`Gene Expression`,
                          min.cells = min.cells,
                          min.features = min.features,
                          assay = "RNA")
  x[["ident"]] <- paste0(ident)
  
  #----- Add complexity measure
  x$log10GenesPerUMI <- log10(x$nFeature_RNA) / log10(x$nCount_RNA)
  return(x)
}