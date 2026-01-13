# Keju
<!-- badges: start -->
<!-- badges: end -->

**keju** is an R package for statistical analysis of Massively Parallel Reporter Assay (MPRA) data, and outputs transcription rate estimates and differential activity estimates with batch-specific uncertainty quantification.

## Installation

### R package installation

**keju** relies on [cmdstanr](https://mc-stan.org/cmdstanr/) and [basilisk](https://www.bioconductor.org/packages/release/bioc/html/basilisk.html). To install cmdstanr, run

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

# use cmdstanr to install CmdStan, this requires a working C++ toolchain and compiler
library(cmdstanr)
install_cmdstan(cores = 2)
```

The compiler requirements can be seen at
[stan-dev](https://github.com/stan-dev/stan/wiki/Coding-Style-and-Idioms#supported-cpp-versions-and-compilers).
If you run into issues with installation, please ensure your gcc version
is \> 5. Additionally, if you run into issues with the TBB library, try
downgrading CmdStan to version 2.33.1. 

**keju** performs some sparse matrix processing steps that require python. As a result, the R API is a thin wrapper around these python calls, which require python==3.13.3, numpy==2.2.5, pandas==2.2.3, and [formulaic](https://matthewwardrop.github.io/formulaic/latest/)==1.1.1. We use basilisk to install a new python environment with these packages. See [basilisk](https://www.bioconductor.org/packages/release/bioc/html/basilisk.html), specifically `basilisk::configureBasiliskEnv()`, for more information and for details about installation location. To install basilisk, run
``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("basilisk")
```

After installing cmdstanr and basilisk, install **keju** by running 

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("pimentellab/keju")
library(keju)
```

## Using keju

An introductory vignette can be found
[here](TODO). If you run into problems, please submit and issue on github or email <albertsxue@gmail.com>.
