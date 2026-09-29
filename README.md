
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FARL: Use Factor-Augmented Regularized Latent Regression for Large-Scale Assessment

<!-- badges: start -->

<!-- badges: end -->

<div class="farl-home-logo-wrap">

<img class="farl-home-logo" src="man/figures/logo.png" width="340" alt="FARL logo" />

</div>

FARL implements factor-augmented regularized latent regression (FARLR)
for large-scale assessments. FARLR addresses high-dimensional, strongly
dependent background variables by decomposing their variation into
common latent factors and idiosyncratic components. Both components
enter the latent regression, while sparse regularization selects a
parsimonious set of relevant covariates. This structure stabilizes
estimation, preserves the interpretability of the selected background
variables, and supports plausible-value generation for group-level
inference.

- `Farlr_mml()`: Provides the unified interface for fitting FARLR
  models. Use `method = "FARLR_EMM"` for importance-sampling-based
  expectation-maximization-maximization estimation, or
  `method = "FARLR_Debias"` to correct shrinkage bias in the regularized
  coefficients.

- `Farlr_drawPVs()`: Draws plausible values from the posterior
  latent-trait distributions implied by a fitted FARLR model.

- `Dire_mml()`: Fits the PCA-based DIRE latent regression used as a
  comparison approach for high-dimensional background variables.

- `Dire_drawPVs()`: Draws plausible values from fitted DIRE models.

<br clear="both" />

## Installation

To install this package from source:

1)  **Windows** users may need to install
    [Rtools](https://CRAN.R-project.org/bin/windows/Rtools/) and select
    the option to add Rtools to the system path. **macOS** users may
    need to install the Xcode command-line tools by running
    `sudo xcode-select --install` in Terminal and then install the [GNU
    Fortran compiler](https://mac.r-project.org/tools/). Most **Linux**
    distributions already include suitable compilers, or provide them
    through the system package manager.

2)  Install `devtools` if necessary, and install FARL from GitHub:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("MAP-LAB-UW/FARL", build_vignettes = TRUE)
torch::install_torch()
```

## Tutorial

After installing the `FARL` package, open its tutorial by running:

``` r
vignette("FARL")
```
