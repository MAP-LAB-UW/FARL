# FARLR Marginal Maximum Likelihood Estimation

Fits a Factor-Adjusted Regularized Latent Regression (FARLR) model using
a marginal maximum likelihood (MML) framework. This function serves as a
unified interface that supports multiple estimation strategies,
including `"FARLR_EMM"` and `"FARLR_Debias"`. Latent traits are
integrated out using Monte Carlo approximation, and regression
parameters are estimated under regularization. Depending on the
specified method, the algorithm either employs an EM–M–type iterative
scheme or a post-selection debiasing procedure.

## Usage

``` r
Farlr_mml(
  X,
  Y,
  parTab,
  n_sam = 15,
  method = "FARLR_EMM",
  lambda = seq(0.1, 0.5, by = 0.1),
  delta.criteria = 0.001,
  iter.max = 500,
  window.size = 50,
  verbose = TRUE
)
```

## Arguments

- X:

  Matrix. Covariate design matrix for the latent regression model,
  typically of dimension `n x p`.

- Y:

  Matrix. Observed item response matrix of dimension `n x J`.

- parTab:

  Data frame. Item parameter table containing at least the columns
  `slope`, `difficulty`, and `guessin`, used to define the item response
  model.

- n_sam:

  Integer. Number of Monte Carlo samples per individual used to
  approximate integrals over latent traits. Defaults to `5`.

- method:

  Character string. Estimation method to be used. Supported values
  include:

  `"FARLR_EMM"`

  :   Iterative FARLR estimation based on an EM–M–type updating scheme
      with regularization.

  `"FARLR_Debias"`

  :   FARLR estimation with regularized variable selection followed by a
      post-selection debiased refit.

  Defaults to `"FARLR_EMM"`.

- lambda:

  Numeric vector. Candidate regularization parameters used for penalized
  regression. Defaults to `seq(0.1, 0.5, by = 0.1)`.

- delta.criteria:

  Numeric. Convergence tolerance for iterative updates. Defaults to
  `1e-3`.

- iter.max:

  Integer. Maximum number of iterations allowed for each value of
  `lambda`. Defaults to `200`.

- window.size:

  Integer. Window size for sliding-window averaging of coefficient
  updates used to stabilize iterative estimation. Defaults to `50`.

- verbose:

  Logical. If `TRUE`, progress messages and iteration status are
  displayed during model fitting. Defaults to `TRUE`.

## Value

A list containing estimation results. The exact contents depend on the
selected `method`, but typically include:

- `coefficients`:

  Estimated regression coefficients.

- `sigma`:

  Estimated residual standard deviation.

- `LogLik`:

  Value of the objective function evaluated at the selected
  regularization parameter.

- `lambda`:

  Selected regularization parameter.

- `Convergence`:

  Indicator of convergence status.

## Details

The function marginalizes over latent variables using Monte Carlo
integration and estimates regression parameters under regularization.
When `method = "FARLR_EMM"`, parameters are updated iteratively using an
EM–M–style procedure. When `method = "FARLR_Debias"`, a penalized
estimator is first used for variable selection, followed by a debiased
refit on the selected active set. The regularization parameter is
selected by minimizing a BIC-type criterion over the supplied `lambda`
grid.

## See also

[`farlr_debias`](https://yijunchenguw.github.io/FARL/reference/farlr_debias.md),
`glmnet`,
[`simdata`](https://philchalmers.github.io/mirt/reference/simdata.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit_emm <- farlr_mml(
  X = X,
  Y = Y,
  parTab = parTab,
  method = "FARLR_EMM",
  verbose = TRUE
)

fit_debias <- farlr_mml(
  X = X,
  Y = Y,
  parTab = parTab,
  method = "FARLR_Debias",
  verbose = TRUE
)
} # }
```
