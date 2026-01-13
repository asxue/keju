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
An introductory vignette can be found [here](TODO). If you run into problems, please submit and issue on github or email <albertsxue@gmail.com>.

### Choosing the correct keju model
**keju** is not a singular model, but a suite of models for different use cases that fit within each other, a little like Russian nesting dolls. As a result, choosing the correct model for your use case can be confusing. 

The most general model, and our default recommended starting point, is `no_motif`. `no_motif` makes no assumptions on correlations between enhancers or architectures, and shrinks transcription rate estimates and effect size estimates towards a generic standard normal prior. Use `no_motif` if you don't have a concrete structure among your tested architectures, which many users will not. Even if you do have some kind of structure, `no_motif` will give you viable estimates, just without some motif-level bells and whistles.

In contrast, some users can use `motif_shrinkage` if they have an exploitable structure in their tested architectures. An example of motif-level structure in the architectures is provided below from the [Zahm et al.](TODO) data, where several architectures are all actually testing the same Rarb transcription factor binding motif. In this case, `motif_shrinkage` shrinks estimates towards a shared motif-level mean, and will also provide *motif-level* estimates for transcription rate and effect sizes, not just architecture-level estimates. For example, in the example below `no_motif` and `motif_shrinkage` will both provide 18 transcription rate estimates and 18 effect size estimates, one for each architecture. However, `motif_shrinkage` will also provide a transcription rate estimate and effect size estimate for the motif itself, and the architecture-level estimates will be regularized to these motif-level estimates. Use `motif_shrinkage` if you have some kind of structure among your architectures. If you do not know if your structure qualifies, we recommend just using `no_motif`.

<p>
<img src="man/figures/motif.png" width="300">
</p>

Finally, the most specialized model is `covariate_motif_slope_intercept`, also the full **keju** model. `covariate_motif_slope_intercept` requires motif-level structure, and is used in the case where the covariates provided in the model have interesting effects on transcription rate that we want to disentangle. As an example, in the [Zahm et al.](TODO) data, choice of minimal promoter strongly affects transcription rate estimates of each architecture (see Figure 4 in our paper). Use `covariate_motif_slope_intercept` if you want covariate-level slope and intercept estimates on transcription rate.

If you have questions, the [vignette](TODO) may be helpful.