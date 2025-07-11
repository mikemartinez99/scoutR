# Single-cell Omics Utilities in R (scoutR)

![status](https://img.shields.io/badge/status-active--development-orange)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![R Version](https://img.shields.io/badge/R-4.4.3-blue)
![version](https://badgen.net/badge/version/0.0.1/blue)

A collection of single-cell omics utility functions to aid in data-wrangling, preprocessing, quality-control and plotting automation to streamline and unify analyses. 

Currently supports scRNA, scATAC, and 10X Multiome. 

<img src="/img/SCoutR_HexLogo.png" width="350px" height="400px" />


## Installation

To intall the scoutR package through R, run the following command:

```r
library(devtools)
devtools::install_github("https://github.com/mikemartinez99/scoutR", force = TRUE)
library(scoutR)

```

## General Helper Functions

|Function|Purpose|
|--------|-------|
|`read_10x_GEX`|Reads in GEX data from 10x into a `Seurat` object, filters for user-specified min.cells and min.features|
|`read_10x_ATAC`|Reads in chromatin-accessibility data into a `Signac` and `Seurat` object, filters for user-specified min.cells and min.features, as well as user-specified blacklist regions, plots density plots of QC metrics|
|`mergeMetadata`|Merge metadata from multiple `Seurat` objects into a single dataframe for plotting purposes|

## Quality Control Visualizations

|Function|Purpose|
|--------|-------|
|`plotOutliers`|Plot violin functions of metadata categories, highlights outlier distributions in red|
|`qcViolin`|Plot violin plots of metadata categories split by `orig.ident`|
|`qcThreshold`|Plot violin or ridge-plots plots of metadata categories split by `orig.ident` with a line to demarcate potential filtering thresholds|
|`qcRidges`|Plot ridge plots of metadata categories split by `orig.ident`|
|`qcScatter`|Plot scatter plot of two metadata cateogories colored by `orig.ident` with distribution annotation and sample facet features`|

## Thresholding and Outlier Detection

|Function|Purpose|
|--------|-------|
|`markOutliers`|Wrapper function for `scuttle::isOutlier` to operate on multiple metadata columns at a time|
|`trackCells`|Diagnostic tool: given a named list of thresholds to test, tracks how many cells you would lose per sample (taking double-dipping cells into account) and outputs as a data frame (no actual filtering occurs here)|
|`calcMAD`|Calculates the median-average deviation for a metadata variable of interest at either the upper or lower bound to help determine quantitative thresholds|

## Demo

For a brief demo of functionalities, see [The scoutR Vignette](https://github.com/mikemartinez99/scoutR/blob/main/vignettes/scoutR.md)

## Change Log

**Version 0.0.1** 
- Added `read_10x_GEX` `read_10x_ATAC` `mergeMetadata` `plotOutliers` `qcViolin` `qcThreshold` `qcRidges` `qcScatter` `markOutliers` `trackCells` and `calcMAD`
- Added basic scoutR vignette (to be expanded upon at a later-date)

## Contact

Michael Martinez M.S.

f007qps@dartmouth.edu
