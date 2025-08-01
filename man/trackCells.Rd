% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trackCells.R
\name{trackCells}
\alias{trackCells}
\title{trackCells}
\usage{
trackCells(sampleList, thresholds)
}
\arguments{
\item{sampleList}{A named list of \code{Seurat} objects with valid meta.data slots}

\item{thresholds}{A named list of expressions to pass to \code{seurat::subset}
(the names will serve as column names in the output dataframe)}
}
\value{
A dataframe
}
\description{
Generate a data-frame showing how many cells you lose for each
specific parameter threshold. This function takes a named list of \code{Seurat}
objects and a named list of expressions to pass to \code{seurat::subset}.
This function DOES NOT directly apply the filtering, It is meant to be used a diagnostic.
This function takes into account barcodes that are "double-dippers" across different
metric (i.e., a barcode that is filtered independently across multiple metrics.)
}
\examples{
# Generate individual plots as well as a merged plot
atacList <- list("obj1", "obj2", "obj3")
thresholds <- list(
  pct_reads_in_peaks = "pct_reads_in_peaks > 50",
  blacklist_fraction = "blacklist_fraction < 0.03",
  TSS_enrichment_low = "TSS.enrichment > 2",
  TSS_enrichment_high = "TSS.enrichment < 15",
  nucleosome_signal = "nucleosome_signal < 4")
# Run the function
cellCounts <- trackCells(
  atacList, 
  thresholds)

}
