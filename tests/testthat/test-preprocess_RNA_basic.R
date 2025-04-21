#----- Load libraries
library(testthat)
library(Seurat)

#------------------------------------------------------------------------------#
# GENERATE SOME DUMMY DATA
#------------------------------------------------------------------------------#

#----- Simulate a small Seurat object
counts <- matrix(
  data = rpois(200, lambda = 5),
  nrow = 20,
  ncol = 10,
  dimnames = list(paste0("gene", 1:20), paste0("cell", 1:10))
)

#----- Add some mitochondrial genes
rownames(counts)[1:3] <- paste0("mt-", 1:3)

#----- Create dummy seurat object
seurat_obj <- CreateSeuratObject(counts = counts)

#------------------------------------------------------------------------------#
# UNIT TESTS!ß
#------------------------------------------------------------------------------#

#----- Test that filtering works
test_that("Cells are filtered based on nCount_RNA", {
  result <- preprocess_RNA_basic(seurat_obj, ident = "test", nCount_thresh = 10, mito_gene_prefix = "^mt-")
  expect_true(all(result$nCount_RNA >= 10))
})

#----- Test that default threshold of 300 is applied for nCount_RNA
test_that("Default threshold of 300 is applied", {
  # artificially inflate some nCount_RNA values
  seurat_obj$nCount_RNA <- c(rep(100, 5), rep(400, 5))
  result <- preprocess_RNA_basic(seurat_obj, ident = "test", mito_gene_prefix = "^mt-")
  expect_true(all(result$nCount_RNA >= 300))
})

#----- Test that metadata columns are added
test_that("percentMT and log10GenesPerUMI are added", {
  result <- preprocess_RNA_basic(seurat_obj, ident = "test", nCount_thresh = 50, mito_gene_prefix = "^mt-")
  expect_true("percentMT" %in% colnames(result[[]]))
  expect_true("log10GenesPerUMI" %in% colnames(result[[]]))
})

#----- Test that the output is a valid Seurat object
test_that("preprocess_RNA_basic returns a Seurat object", {
  result <- preprocess_RNA_basic(seurat_obj, 
                                 ident = "TestSample", 
                                 nCount_thresh = 50, 
                                 mito_gene_prefix = "^mt-")
  
  expect_s4_class(result, "Seurat")
})

#----- Test that error messages will work as expected
test_that("errors messages work", {
  expect_error(
    preprocess_RNA_basic(seurat_obj, ident = "Test1", nCount_thresh = 300, mito_gene_prefix = "^mt-"),
    paste0("No cells found")
  )
})

#----- Test overall function
result <- preprocess_RNA_basic(seurat_obj, ident = "Test1", nCount_thresh = 100, mito_gene_prefix = "^mt-")
resultMeta <- result@meta.data
