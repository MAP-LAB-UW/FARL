# Simulated Dataset: Multidimensional FARLR Item Responses

A simulated multidimensional item-response dataset for illustrating and
testing multidimensional factor-augmented regularized latent regression
(MFARLR) and multidimensional DIRE. The dataset contains five correlated
latent dimensions, high-dimensional covariates, and simple-structure
dichotomous item responses.

## Usage

``` r
data(sim_m1)
```

## Format

A list with the following components:

- X:

  Numeric covariate matrix with one row per respondent.

- Y:

  Dichotomous item-response matrix with one row per respondent and one
  column per item.

- theta:

  Matrix of true multidimensional latent-trait values.

- a:

  Item discrimination matrix. Rows represent items and columns represent
  latent dimensions.

- b:

  Vector of true item difficulty parameters.

- d:

  Vector of item intercept parameters used for response simulation.

- parTab:

  Item parameter table used by the model-fitting functions.

- item_index:

  List identifying the items assigned to each latent dimension.

- main:

  Integer indices of the covariates treated as main variables.

- main_vars:

  Character names of the main covariates.

- X_discrete:

  Matrix containing the discretized main covariates.

- beta:

  Matrix of true common-factor regression coefficients.

- veta:

  Matrix of true idiosyncratic regression coefficients.

- active_index:

  Indices of the truly active covariates for each latent dimension.

- Sigma_error:

  Residual covariance matrix used before latent-trait standardization.

- true_prior_mean:

  Matrix containing the true conditional latent means on the
  standardized latent-trait scale.

- true_sigma:

  True conditional latent covariance matrix on the standardized
  latent-trait scale.

- theta_center:

  Centering constants used to standardize the latent dimensions.

- theta_scale:

  Scaling constants used to standardize the latent dimensions.

## Details

The covariates are generated from a two-factor model with additional
idiosyncratic noise. Each latent dimension depends on the common factors
and a sparse set of covariates. Correlated multivariate residual errors
are added before the latent dimensions are centered and standardized.
Items follow a simple loading structure so that each item measures one
latent dimension.

## Author

Yijun Cheng <chengxb@uw.edu>
