# FARL: Factor-Augmented Regularized Latent Regression for Large-Scale Assessment

FARL is developed to support large-scale assessment (LSA) analyses, in
which latent regression models are used to integrate students’
background information and to generate plausible values (PVs) for
secondary analysis. In LSA settings, background variables are often
high-dimensional and highly correlated, posing substantial challenges
for traditional latent regression approaches.

## Details

FARL implements factor-augmented regularized latent regression (FARLR),
an innovative framework that jointly models common factors and
idiosyncratic components. Regularization is employed to select relevant
idiosyncratic predictors, yielding interpretable regression results
while maintaining congeniality and estimation stability.

## Latent regression models

- [`Farlr_mml`](https://yijunchenguw.github.io/FARL/reference/Farlr_mml.md)
  fits the factor-augmented regularized latent regression model via
  marginal maximum likelihood

  - Method 1: `FARLR_EMM`

  - Method 2: `FARLR_Debias`

&nbsp;

- [`Dire_mml`](https://yijunchenguw.github.io/FARL/reference/Dire_mml.md)
  fits the DIRE latent regression model

## Plausible value generation

- [`Farlr_drawPVs`](https://yijunchenguw.github.io/FARL/reference/Farlr_drawPVs.md)
  generates plausible values under the FARLR framework

&nbsp;

- [`Dire_drawPVs`](https://yijunchenguw.github.io/FARL/reference/Dire_drawPVs.md)
  generates plausible values under the Dire framework

## Simulation and example data

- Built-in example datasets illustrating FARLR estimation and PV
  generation

- Utility functions for simulation studies in large-scale assessment
  settings

## See also

Useful links:

- <https://yijunchenguw.github.io/FARL/>

- <https://github.com/yijunchenguw/FARL>

## Author

**Maintainer**: Yijun Cheng <chengxb@uw.edu>
([ORCID](https://orcid.org/0000-0002-0671-9193))
