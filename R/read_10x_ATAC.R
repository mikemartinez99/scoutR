#' @title read_10x_ATAC
#' 
#' @description
#' When multiple samples are present, `read_10x_ATAC` can read each sample's 
#' filtered feature barcode matrix h5 file and generate a `ChromatinAssay` and 
#' `Seurat` object for each. Each of these objects is then combined into a list 
#' of objects. This function can be easily parallelized using `purrr::map`. 
#' This function requires that you have your organism annotations as a `GRanges`
#' object and blacklist region file in bed format. By default, it will output 
#' density plots of nCount peaks vs TSS enrichment.
#' 
#' @param sample A sample name (must be same name as your cellranger outputs to dynamically construct file paths)
#' @param files_dir A directory path ending in "/" where raw 10x data lives
#' @param annotations A GRanges object with organism annotations 
#' @param min.cells Minimum number of cells a feature needs to be present in to be retained
#' @param min.features Minimum number of features needed in a cell to be retained
#' @param blacklist A blacklist region file in bed format 
#' @param figDir A directory path ending in "/" where density plots will be output to
#'
#' @returns A list of seurat objects (one per sample)
#' 
#' @examples # Path to sample directory and output directory
#' wd <- "/my/path/to/data/"
#' outDir <- "/my/path/to/data/outputs/"
#' # Vector of sample IDs
#' sampleIDs <- c("S1", "S2", "S3")
#' # Vector of idents
#' idents <- c("A", "B", "C")
#' # Read in annotations and blacklist bed
#' anno <- readRDS("path/to/annotations.Rds")
#' blacklist <- read.csv("/path/to/blacklist.bed", sep = "\t")
#' # Call the function
#' seurat_list <- map(samples, 
#'   ~ run_atac_qc_filtered(
#'     .x, 
#'     wd, 
#'     anno, 
#'     min.cells = 10, 
#'     min.features = 200, 
#'     blacklist, 
#'     outDir)
#'     )
#' 
#' @importFrom purrr map 
#' @import Seurat
#' @import Signac
#' @import GenomicRanges
#' 
#' @export


read_10x_ATAC <- function(sample, files_dir, annotations, min.cells, min.features, blacklist, figDir){
  message(paste0("running sample ", sample))
  #----- CONSTRUCT FILE PATH TO FEATURE MATRIX
  counts <- Read10X_h5(filename = paste0(files_dir, sample, "/outs/filtered_feature_bc_matrix.h5"))
  #----- ADD ANNOTATIONS
  annotations <- annotations
  #-----CREATE CHROMATIN ASSAY 
  chrom_assay <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = paste0(files_dir, sample, '/outs/atac_fragments.tsv.gz'),
    annotation = annotations,
    min.cells = min.cells,
    min.features = min.features)
  #----- READ METADATA AND SUBSET DATA ACCORDINGLY 
  metadata <- read.csv(
    file = paste0(files_dir, sample, "/outs/per_barcode_metrics.csv"),
    header = TRUE)
  rownames(metadata) <- metadata[,1]
  any(is.na(rownames(metadata)))
  any(is.na(colnames(chrom_assay)))
  message(paste0("Number of columns in chrom assay: ", ncol(chrom_assay)))
  message(paste0("Number of rows in metadata: ", nrow(metadata)))
  metadata_sub <- metadata[rownames(metadata) %in% colnames(chrom_assay),]
  #----- CREATE SEURAT OBJECT AND ADD ANNOTATIONS 
  seurat <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks", 
    meta.data = metadata_sub)
  seurat$sample <- paste0(sample)
  Annotation(seurat) <- annotations
  #----- CALCULATE QC METRICS 
  seurat$pct_reads_in_peaks <- seurat$atac_peak_region_fragments / seurat$atac_fragments * 100
  overlaps <- findOverlaps(query = seurat[['peaks']], 
                           subject = blacklist)
  hit.regions <- queryHits(x = overlaps)
  data.matrix <- GetAssayData(object = seurat, 
                              assay = 'peaks', 
                              slot = "counts")[hit.regions, , drop = FALSE]
  seurat$blacklist_region_fragments <- colSums(data.matrix)
  seurat$blacklist_fraction <- seurat$blacklist_region_fragments / seurat$atac_peak_region_fragments
  seurat <- NucleosomeSignal(object = seurat)
  seurat$nucleosome_group <- ifelse(seurat$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  seurat <- TSSEnrichment(object = seurat, fast = FALSE)
  seurat$high.tss <- ifelse(seurat$TSS.enrichment > 2, 'High', 'Low')
  #----- PLOT DENSITY 
  densityPlot <- DensityScatter(seurat,
                                x = "nCount_peaks",
                                y = "TSS.enrichment",
                                log_x = TRUE,
                                quantiles = TRUE)
  fileName <- paste0(sample, "_Density_Plot.png")
  filePath <- paste0(figDir, fileName)
  ggsave(filePath, densityPlot, width = 12, height = 10)
  #----- RETURN SEURAT OBJECT
  return(seurat)
}
