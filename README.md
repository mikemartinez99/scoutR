# Single-cell Omics Utilities in R (scoutR)

![status](https://img.shields.io/badge/status-in--development-orange)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![R Version](https://img.shields.io/badge/R-4.4.3-blue)

A collection of utility and plotting functions to streamline the processing of single-cell omics data (supports scRNA, scATAC, 10X multiome, and Visium)

<img src="/img/SCoutR_HexLogo.png" width="350px" height="400px" />

# 🛠️ Documentation under construction...stay tuned!


## Installation

To intall the scoutR package through R, run the following command:

```r
library(devtools)
devtools::install_github("https://github.com/mikemartinez99/scoutR", force = TRUE)
library(scoutR)

```

## Functionality

- Add options to plotting functions to color/group by another metadata feature besides orig.ident, but keep orig.ident as the default

|Function|Purpose|
|--------|-------|
|`read_10x_GEX`|Reads in GEX data from 10x into a `Seurat` object, filters for user-specified min.cells and min.features|
|`read_10x_ATAC`|Reads in chromatin-accessibility data into a `Signac` and `Seurat` object, filters for user-specified min.cells and min.features, as well as user-specified blacklist regions, plots density plots of QC metrics|
|`mergeMetadata`|Merge metadata from multiple `Seurat` objects into a single dataframe for plotting purposes|
|`markOutliers`|Wrapper function for `scuttle::isOutlier` to operate on multiple metadata columns at a time|
|`plotOutliers`|Plot violin functions of metadata categories, highlights outlier distributions in red|
|`qcViolin`|Plot violin plots of metadata categories split by `orig.ident`|
|`qcThreshold`|Plot violin plots of metadata categories split by `orig.ident` with a horizontal line to indicate thresholds|
|`qcRidges`|Plot ridge plots of metadata categories split by `orig.ident`|
|`qcScatter`|Plot scatter plot of two metadata cateogories colored by `orig.ident` with distribution annotation and sample facet features`|
|`trackCells`|Diagnostic tool: given a named list of thresholds to test, tracks how many cells you would lose per sample and outputs as a data frame (no actual filtering occurs here)|


## Contact

Michael Martinez M.S.

f007qps@dartmouth.edu
