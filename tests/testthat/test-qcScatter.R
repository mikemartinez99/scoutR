

test_that("qcScatter returns ggplot object and expected plot", {
  
  
  # Dummy metadata
  metadata <- data.frame(
    orig.ident = rep(c("SampleA", "SampleB", "SampleC"), each = 10),
    nFeature_peaks = c(rnorm(10, mean = 1000), 
                       rnorm(10, mean = 2000), 
                       rnorm(10, mean = 2000)),
    nCount_peaks = c(rnorm(10, mean = 1000), 
                     rnorm(10, mean = 2000), 
                     rnorm(10, mean = 2000)))
  
  # Dummy colors and temp fig dir
  test_colors <- c("SampleA" = "#00AFBB", "SampleB" = "#E7B800", "SampleC" = "red")  
  figDir <- "/users/mike/Desktop/"
  
  # Call the function
  p <- qcScatter(metadata, Xvar = "nFeature_peaks", Yvar = "nCount_peaks", logTransformX = FALSE, logTransformY = FALSE, test_colors, figDir, width = 4, height = 3)
  
  # Tests
  expect_s3_class(p, "ggplot")
  
  # Check file saved
  expected_file <- file.path(figDir, "nFeature_peaks_vs_nCount_peaks.png")
  expect_true(file.exists(expected_file))
})
