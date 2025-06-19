test_that("qcViolin returns a ggplot object and saves file", {
  
  
  # Dummy metadata
  metadata <- data.frame(
    orig.ident = rep(c("SampleA", "SampleB", "SampleC"), each = 10),
    nFeature_peaks = c(rnorm(10, mean = 1000), 
                       rnorm(10, mean = 2000), 
                       rnorm(10, mean = 2000)))
  
  # Dummy colors and temp fig dir
  test_colors <- c("SampleA" = "#00AFBB", "SampleB" = "#E7B800", "SampleC" = "red")  
  figDir <- "/users/mike/Desktop/"
  
  # Call the function
  p <- qcViolin(metadata, "nFeature_peaks", logTransform = FALSE, test_colors, figDir, width = 4, height = 3)
  
  # Tests
  expect_s3_class(p, "ggplot")
  
  # Check file saved
  expected_file <- file.path(figDir, "Log_nFeature_peaks_QC_Violins.png")
  expect_true(file.exists(expected_file))
})


test_that("qcViolin log transformation works", {
  
  
  # Dummy metadata
  metadata <- data.frame(
    orig.ident = rep(c("SampleA", "SampleB", "SampleC"), each = 10),
    nFeature_peaks = c(rnorm(10, mean = 1000), 
                       rnorm(10, mean = 2000), 
                       rnorm(10, mean = 2000)))
  
  # Dummy colors and temp fig dir
  test_colors <- c("SampleA" = "#00AFBB", "SampleB" = "#E7B800", "SampleC" = "red")  
  figDir <- "/users/mike/Desktop/"
  
  # Call the function
  p <- qcViolin(metadata, "nFeature_peaks", logTransform = TRUE, test_colors, figDir, width = 4, height = 3)
  
  # Tests
  expect_s3_class(p, "ggplot")
  
  # Check file saved
  expected_file <- file.path(figDir, "Log_nFeature_peaks_QC_Violins.png")
  expect_true(file.exists(expected_file))
})

