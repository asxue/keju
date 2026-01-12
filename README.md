# Keju
<!-- badges: start -->
[![R-CMD-check](https://github.com/asxue/keju/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/asxue/keju/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Keju** is an R package for scoring FACS-based DMS experiments with
uncertainty quantification. It takes in a negative control group
(usually synonymous variants) and scores each variant relative to the
negative control group.

## Installation

### R package installation

**Keju** relies on [cmdstanr](https://mc-stan.org/cmdstanr/), which should
be properly installed first.

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


After installing cmdstanr, install **keju** by running 

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("pimentellab/keju")
library(keju)
```
