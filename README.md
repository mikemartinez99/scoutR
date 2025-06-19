# Single-cell Omics Utilities in R (scoutR)

![status](https://img.shields.io/badge/status-in--development-orange)
![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![R Version](https://img.shields.io/badge/R-4.4.3-blue)

Package containing helpful utilities and plotting functions for single cell omics analysis (scRNA, scATAC, 10X multiome)

<img src="/img/SCoutR_HexLogo.png" width="350px" height="400px" />

# Table of Contents
- [Installation](#installation)
- [Functions](#functions)
- [Vignettes](#vignettes)
- [Contact](#contact)

## Installation

To intall the scoutR package through R, run the following command:

```r
library(devtools)
devtools::install_github("mmartinez99/scoutR/")

library(scoutR)

```

## Functions
|Function|Purpose|
|--------|-------|
|`read_10x_multi`|When multiple 10x outputs are present in a directory, read in the filtered feature barcode matrices and generate seurat object for each sample. Output is saved as a list of seurat objects with unique names|
|`preprocess_RNA_basic`|Filter seurat objects based on nCount_RNA threshold and add additional metadata (percent mitochondria and log10 genes per umi)|

## Vignettes

## Contact

Michael Martinez M.S.

f007qps@dartmouth.edu
