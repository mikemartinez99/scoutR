#' @title mergeMeta
#' 
#' @description
#' Quick way to merge metadata actoss multiple Seurat objects into one master dataframe for plotting purposes
#' 
#' @param x list of Seurat objects (could be named or un-named, but a named list is recommended)
#'
#' @returns A dataframe of merged metadata
#' @export
#' @import Seurat
#'
#' @examples # list of seurat objects with naems
#' @examples objects <- c("obj1", "obj2", "obj3", "obj4")
#' @examples # Vector of sample IDs
#' @examples sampleIDs <- c("S1", "S2", "S3", "S4")
#' @examples names(objects) <- sampleIDs
#' @examples metadata <- mergeMetadata(objects)

mergeMetadata <- function(x) {
  
  #----- Check that sample names match length of list
  if (!is.list(seurat_list)) {
    stop("Input must be a list of Seurat objects.")
  } 
  
  #----- Assign sample names if not provided
  sample_names <- names(seurat_list)
  if (is.null(sample_names) || any(sample_names == "")) {
    sample_names <- paste0("sample", seq_along(seurat_list))
  }
  
  #----- Merge metadata
  merged_metadata <- do.call(rbind, lapply(seq_along(x), function(i) {
    meta <- x[[i]]@meta.data
    meta$orig.ident <- sample_names[i]
    meta$cell_id <- rownames(meta)
    return(meta)
  }))
  
  #----- Return metadata
  return(merged_metadata)
}