test_that("calcMAD returns numeric value for upper and lower bounds", {
  # Simulate a basic Seurat object
  mat <- matrix(rpois(100, lambda = 10), nrow = 10)
  seurat_obj <- Seurat::CreateSeuratObject(counts = mat)
  seurat_obj$TSS.enrichment <- c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  
  # Test upper bound
  upper_bound <- calcMAD(
    x = seurat_obj,
    variable = "TSS.enrichment",
    nmads = 2,
    bound = "higher"
  )
  
  expect_type(upper_bound, "double")
  expect_length(upper_bound, 1)
  
  # Test lower bound
  lower_bound <- calcMAD(
    x = seurat_obj,
    variable = "TSS.enrichment",
    nmads = 2,
    bound = "lower"
  )
  
  expect_type(lower_bound, "double")
  expect_length(lower_bound, 1)
  
  # Lower bound should be greater than upper bound in this case
  expect_gt(lower_bound, upper_bound)
})

test_that("calcMAD errors when variable does not exist", {
  mat <- matrix(rpois(100, lambda = 10), nrow = 10)
  seurat_obj <- Seurat::CreateSeuratObject(counts = mat)
  
  expect_error(
    calcMAD(seurat_obj, variable = "nonexistent", nmads = 2),
    "not found in `@meta.data`"
  )
})

test_that("calcMAD errors on non-Seurat input", {
  expect_error(
    calcMAD(x = data.frame(x = 1:10), variable = "x", nmads = 2),
    "`x` must be a Seurat object"
  )
})