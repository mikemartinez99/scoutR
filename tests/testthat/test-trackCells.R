test_that("trackCells returns a data.frame with expected columns", {
  # Simulate input Seurat object(s)
  mat <- matrix(rpois(100, lambda = 10), nrow = 10)
  seurat_obj <- Seurat::CreateSeuratObject(counts = mat)
  seurat_obj$pct_reads_in_peaks <- runif(10, 0, 100)
  seurat_obj$blacklist_fraction <- runif(10, 0, 0.1)
  
  sampleList <- list(S1 = seurat_obj)
  
  thresholds <- list(
    pct_reads_in_peaks = "pct_reads_in_peaks > 50",
    blacklist_fraction = "blacklist_fraction < 0.05"
  )
  
  res <- trackCells(sampleList, thresholds)
  
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("Sample", "Original_Count", "pct_reads_in_peaks", "blacklist_fraction", "Total_Cells_Removed"))
  expect_equal(res$Sample, "S1")
})


test_that("trackCells correctly handles double-dipped cells (not double-counted)", {
  library(Seurat)
  
  # Simulate Seurat object with 5 cells
  mat <- matrix(rpois(50, lambda = 10), nrow = 10)
  seurat_obj <- suppressMessages(CreateSeuratObject(counts = mat))
  
  # Add metadata where 3 cells fail both filters
  seurat_obj$pct_reads_in_peaks <- c(30, 40, 45, 90, 95)  # First 3 < 50
  seurat_obj$blacklist_fraction <- c(0.1, 0.07, 0.06, 0.01, 0.02)  # First 3 > 0.05
  
  sampleList <- list(S1 = seurat_obj)
  
  thresholds <- list(
    pct_reads_in_peaks = "pct_reads_in_peaks > 50",            # Fails for cells 1–3
    blacklist_fraction = "blacklist_fraction < 0.05"           # Fails for cells 1–3
  )
  
  res <- trackCells(sampleList, thresholds)
  
  # Check expected filtered counts per filter
  expect_equal(res$pct_reads_in_peaks, 3)
  expect_equal(res$blacklist_fraction, 3)
  
  # Total cells removed should only be 3, not 6
  expect_equal(res$Total_Cells_Removed, 3)
})
