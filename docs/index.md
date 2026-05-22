# FARL: Use Factor-augmented Regularized Latent Regression for Large Scale Assessment

FARL models integrate latent trait estimation with regularized
regression, enabling scalable modeling of relationships between latent
abilities and high-dimensional covariates.

- FARLR_EMM: Factor-augmented regularized latent regression estimated
  via an iterative EMM-style procedure.

- FARLR_debias: Factor-augmented regularized latent regression with post
  debiasing.

- FARLR_drawPVs: Drawing plausible values from fitted FARLR models.

- Dire_mml: Marginal maximum likelihood estimation of DIRE models.

- Dire_drawPVs: Drawing plausible values from fitted DIRE models.

## Installation

To install this package from source:

1.  **Windows** users may need to install the
    [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/) and include
    the checkbox option of installing Rtools to their path for easier
    command line usage. **Mac** users may need to install Xcode command
    line tools by `sudo xcode-select --install` in the terminal, and
    then install [GNU Fortran
    compiler](https://mac.r-project.org/tools/). Most **Linux**
    distributions should already have up to date compilers (or if not
    they can be installed/updated easily).

2.  Install the `devtools` package (if necessary), and install the
    package from [GitHub](https://github.com/) with

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yijunchenguw/FARL", build_vignettes = TRUE)
torch::install_torch()
```

## Tutorial

After installing the `FARL` package, open its tutorial by running

``` r
vignette("FARL")
```
