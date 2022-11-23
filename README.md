
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mGWASR <img src="man/figures/mgwas_logo.png" align="right" alt="" width="120" />

#### A companion R package for the mGWAS-Explorer web server

<!-- badges: start -->
<!-- badges: end -->

Note - mGWASR is still under development - we cannot guarantee full
functionality.

## Description

**mGWASR** contains the R functions and libraries underlying the
[mGWAS-Explorer](https://www.mgwas.ca/) web server. After installing and
loading the package, users will be able to reproduce the same results
from their local computers using the corresponding R command history
downloaded from mGWAS-Explorer, thereby achieving maximum flexibility
and reproducibility.

## Getting Started

### Step 1. Install package dependencies

To use mGWASR , first install all package dependencies. Ensure that you
are able to download packages from bioconductor. To install package
dependencies, use the pacman R package (for those with \>R 3.5.1). Note
that some of these packages may require additional library dependencies
that need to be installed prior to their own successful installation.

``` r
install.packages("pacman")

library(pacman)

# need to update the depdencies

pacman::p_load(RSQLite, igraph,  BiocParallel,  pryr,  httr,  reshape,  ggplot2,  RJSONIO,  RCurl,  XML,  ggforce,  graphlayouts,  compiler,  dplyr,  RColorBrewer,  Cairo,  plyr,  qs,  rjson,  TwoSampleMR)
```

### Step 2. Install the package

mGWASR is freely available from GitHub. The package documentation,
including the vignettes for each module and user manual is available
within the downloaded R package file. If all package dependencies were
installed, you will be able to install the mGWASR. Due to issues with
Latex, some users may find that they are only able to install mGWASR
without any documentation (i.e. vignettes).

Install the package directly from github using the *devtools* package.
Open R and enter:

``` r
# Step 1: Install devtools
install.packages("devtools")
library(devtools)

# Step 2: Install mGWASR WITHOUT documentation
devtools::install_github("xia-lab/mGWASR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"), force = TRUE)

# Step 2: Install mGWASR WITH documentation
devtools::install_github("xia-lab/mGWASR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE, force = TRUE)
```

## Tutorials

For detailed tutorials on how to use mGWASR, please refer to the R
package vignettes. Note, the functions below work only if the R package
vignettes were built.

Within R:

``` r
vignette(package="mGWASR")
```

Within a web-browser:

``` r
browseVignettes("mGWASR")
```

## Citation

mGWASR has been developed by the [Xia Lab](http://xialab.ca/) at McGill
University. The original manuscript (web-based version) can be found
[here](https://www.mdpi.com/2218-1989/12/6/526).

We encourage users to further develop the package to suit their needs.
If you use the R package, please cite us:

Chang, L., Zhou, G., Ou, H., & Xia, J. (2022). mGWAS-Explorer: Linking
SNPs, Genes, Metabolites, and Diseases for Functional Insights.
Metabolites, 12(6), 526. <https://doi.org/10.3390/metabo12060526>

## Bugs or feature requests

To inform us of any bugs or requests, please open a new issue or send an
email to le.chang @ mail.mcgill.ca
