# FARL: Use Factor-Augmented Regularized Latent Regression for Large-Scale Assessment

FARL supports latent regression analyses for large-scale assessments
(LSAs), where student background information is incorporated into the
population model used to generate plausible values (PVs). These
background variables are often high-dimensional and strongly dependent,
which can make conventional latent regression unstable and complicate
variable selection.

## Details

FARL implements factor-augmented regularized latent regression (FARLR).
FARLR decomposes covariate variation into common latent factors and
idiosyncratic components, includes both components in the latent
regression, and uses sparse regularization to select a parsimonious set
of relevant covariates. This structure is designed to stabilize
estimation while retaining interpretable background-variable effects.

## Model estimation

- [`Farlr_mml`](https://map-lab-uw.github.io/FARL/reference/Farlr_mml.md)
  is the unified interface for fitting FARLR models. It supports two
  estimation methods:

  - `method = "FARLR_EMM"`: importance-sampling-based
    expectation-maximization-maximization estimation, with a regularized
    variable-selection step followed by an unpenalized coefficient
    update.

  - `method = "FARLR_Debias"`: regularized estimation using the factor
    and idiosyncratic components, followed by correction of LASSO
    shrinkage bias.

- [`Dire_mml`](https://map-lab-uw.github.io/FARL/reference/Dire_mml.md)
  fits a PCA-based DIRE latent regression model for comparison with the
  FARLR approaches.

- [`MFARLR_mml`](https://map-lab-uw.github.io/FARL/reference/MFARLR_mml.md)
  fits multidimensional FARLR by combining dimension-specific sparse
  covariate screening with joint Gaussian variational
  expectation-maximization estimation.

- [`MDIRE_mml`](https://map-lab-uw.github.io/FARL/reference/MDIRE_mml.md)
  fits the multidimensional PCA-based DIRE comparison model.

## Plausible-value generation

- [`Farlr_drawPVs`](https://map-lab-uw.github.io/FARL/reference/Farlr_drawPVs.md)
  generates plausible values from the posterior latent-trait
  distributions implied by a fitted FARLR model.

- [`Dire_drawPVs`](https://map-lab-uw.github.io/FARL/reference/Dire_drawPVs.md)
  generates plausible values from a fitted DIRE model.

- [`MFARLR_drawPVs`](https://map-lab-uw.github.io/FARL/reference/MFARLR_drawPVs.md)
  generates correlated multidimensional plausible values from a fitted
  MFARLR model.

- [`MDIRE_drawPVs`](https://map-lab-uw.github.io/FARL/reference/MDIRE_drawPVs.md)
  generates plausible values from a fitted multidimensional DIRE model.

## Simulation and example data

- [`sim_a1`](https://map-lab-uw.github.io/FARL/reference/sim_a1.md)
  contains simulated dichotomous 2PL responses, item parameters, and
  high-dimensional background covariates.

- [`sim_a2`](https://map-lab-uw.github.io/FARL/reference/sim_a2.md)
  contains simulated mixed 3PL and GPCM responses, item parameters, and
  high-dimensional background covariates.

- [`sim_m1`](https://map-lab-uw.github.io/FARL/reference/sim_m1.md)
  contains simulated multidimensional item responses, item parameters,
  latent traits, and background covariates.

## Methodological scope

FARLR is intended to improve compatibility between plausible-value
generation and common secondary analyses, but it does not guarantee
congeniality in the strict statistical sense. The package supports both
unidimensional and multidimensional latent proficiency models;
currently, `MFARLR_mml` supports multidimensional binary 2PL items.

## See also

Useful links:

- <https://MAP-LAB-UW.github.io/FARL/>

- <https://github.com/MAP-LAB-UW/FARL>

## Author

**Maintainer**: Yijun Cheng <chengxb@uw.edu>
([ORCID](https://orcid.org/0000-0002-0671-9193))

Authors:

- Chun Wang <wang4066@uw.edu>
  ([ORCID](https://orcid.org/0000-0003-2695-9781))

- Gongjun Xu <gongjun@umich.edu>
  ([ORCID](https://orcid.org/0000-0003-4023-5413))
