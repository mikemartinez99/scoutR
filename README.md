## IN DEVELOPMENT ##

# Single-cell Omics Utilities in R (scoutR)
This package contains utility functions to help navigate and preprocess single-cell omics data (RNA, ATAC, and 10X Multiome). 

**Package creator:** Mike Martinez M.S.

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
