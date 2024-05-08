# smoothclust

[![R build status](https://github.com/lmweber/smoothclust/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/smoothclust/actions)


## Overview

`smoothclust` is a method for segmentation of spatial domains and spatially-aware clustering in spatial transcriptomics data. The method generates spatial domains with smooth boundaries by smoothing gene expression profiles across neighboring spatial locations, followed by unsupervised clustering. Spatial domains consisting of consistent mixtures of cell types may then be further investigated by applying cell type compositional analyses or differential analyses.


## Installation

The `smoothclust` package can be installed from Bioconductor as follows (using R version 4.4 onwards). This is the recommended installation for most users. Additional details are shown on the [Bioconductor](https://bioconductor.org/packages/smoothclust) package landing page.

```
install.packages("BiocManager")
BiocManager::install("smoothclust")
```

The latest development version of the package can also be installed from the [development version of Bioconductor](https://contributions.bioconductor.org/use-devel.html), or alternatively from GitHub as follows.

```
remotes::install_github("lmweber/smoothclust")
```


## Tutorial

For a tutorial and example workflow, see the package vignette.


## Citation

In preparation.
