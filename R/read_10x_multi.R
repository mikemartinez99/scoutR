#' @title read_10x_multi
#' 
#' @description
#' When multiple samples are present, `read_10x_multi` can read each sample's filtered feature barcode matrix h5 file
#' and generate a seurat object for each. Each of these seurat objects is then combined into a list of objects.
#' This function can be easily parallelized using `purrr::map`
#' 
#' @importFrom purrr map 
#' @import Seurat
#'
#' @param sampleDir Path to the sample directory containing multiple subfolders (one per sample)
#' @param sample A vector of sample names that correspond to how subfolders in `sampleDir` are named
#' @param ident A vector of desired ident names for seurat object. Idents must correspond to same order in `sample`
#'
#' @returns A list of seurat objects (one per sample)
#' @export
#'
#' @examples # Path to sample directory
#' @examples wd <- "/my/path/to/data/"
#' @examples # Vector of sample IDs
#' @examples sampleIDs <- c("S1", "S2", "S3")
#' @examples # Vector of idents
#' @examples idents <- c("A", "B", "C")
#' @examples # Call the function
#' @examples seuratList <- purrr::map(sampleIDs, read_10x_multi(wd, sampleIDs, idents))

read_10x_multi <- function(sampleDir, sample, ident) {
  print(sample)
  counts <- Read10X_h5(filename = paste0(sampleDir, sample, "/outs/filtered_feature_bc_matrix.h5"))
  x <- CreateSeuratObject(counts = counts$`Gene Expression`,
                          min.cells = 10,
                          min.features = 100,
                          assay = "RNA")
  x[["ident"]] <- paste0(ident)
  return(x)
}