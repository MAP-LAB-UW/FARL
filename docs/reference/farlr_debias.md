# FARLR Debiased Estimation for Regularized Latent Regression

Fits a Factor-Augmented Regularized Latent Regression (FARLR) model with
a post-selection debiasing procedure. Latent traits are approximated
using a Monte Carlo scheme based on a normal approximation to IRT
posterior estimates (`theta_est_irt.mean`, `theta_est_irt.se`).
Regression coefficients are first obtained via weighted LASSO
regularization using `glmnet`, after which a debiased refit is performed
on the selected active set using weighted least squares. To improve
numerical stability, coefficient updates are smoothed across iterations
using a sliding-window averaging scheme. The regularization parameter
`lambda` is selected by minimizing a BIC-type criterion over
`lambda_all`.

## Usage

``` r
farlr_debias(
  n,
  resp,
  parTab,
  K_hat,
  p,
  lambda_all,
  delta.criteria = 0.001,
  iter.max = 500,
  n_sam = 5,
  window.size = 50,
  theta_est_irt.mean,
  theta_est_irt.se,
  resp_rep = NA,
  Z.em = NA,
  Uupdate,
  hatU,
  Fan,
  main,
  verbose = TRUE
)
```

## Arguments

- n:

  Integer. Number of individuals (sample size).

- resp:

  Matrix. Observed item response matrix of dimension `n x J`.

- parTab:

  Data frame. Item parameter table containing at least the columns
  `slope`, `difficulty`, and `guessin`.

- K_hat:

  Integer. Number of latent components included in the regression design
  matrix `Z.em`.

- p:

  Integer. Number of observed covariates included in `Z.em`.

- lambda_all:

  Numeric vector. Candidate regularization parameters supplied to
  `glmnet`.

- delta.criteria:

  Numeric. Convergence tolerance for iterative updates. Defaults to
  `1e-3`.

- iter.max:

  Integer. Maximum number of iterations for each value of `lambda`.
  Defaults to `500`.

- n_sam:

  Integer. Number of Monte Carlo samples per individual used to
  approximate latent trait uncertainty. Defaults to `5`.

- window.size:

  Integer. Window size for sliding-window averaging of regression
  coefficient updates. Defaults to `50`.

- theta_est_irt.mean:

  Numeric vector of length `n`. Posterior mean estimates of the latent
  trait obtained from an IRT model.

- theta_est_irt.se:

  Numeric vector of length `n`. Posterior standard error estimates of
  the latent trait obtained from an IRT model.

- resp_rep:

  Matrix. Replicated or expanded response matrix used for Monte Carlo
  integration (see `q_num_NA`).

- Z.em:

  Matrix. Regression design matrix, typically of dimension
  `n x (K_hat + p)`.

- Uupdate:

  Internal object. Passed to internal update routines controlling
  coefficient smoothing.

- hatU:

  Internal object. Passed to internal routines used in the debiasing
  refit.

- Fan:

  Internal object. Passed to internal routines for adaptive weighting or
  regularization.

- main:

  Integer vector. Indices of predictors to be treated as unpenalized in
  `glmnet` via `penalty.factor = 0`. Remaining predictors are penalized
  unless they correspond to the first `K_hat` latent components.

- verbose:

  Logical. If `TRUE`, progress messages and a progress indicator are
  displayed. Defaults to `TRUE`.

## Value

A list containing:

- `coefficients`:

  Debiased regression coefficient estimates corresponding to the
  selected `lambda` (length `K_hat + p`).

- `sigma`:

  Estimated residual standard deviation.

- `LogLik`:

  Value of the BIC-type objective function at the selected `lambda`.
  Returned as `LogLik` for compatibility.

- `minBIC`:

  Index of `lambda_all` that minimizes the BIC-type criterion.

- `Convergence`:

  Character string indicating convergence status of the iterative
  procedure.

## Details

For each candidate value in `lambda_all`, coefficient estimates are
updated iteratively until convergence or until `iter.max` iterations are
reached. Convergence is assessed using the maximum absolute change in
regression coefficients and the residual scale parameter: \$\$ \delta =
\max\left( \lvert \sigma^{(t)} - \sigma^{(t-1)} \rvert, \max_j \lvert
\beta_j^{(t)} - \beta_j^{(t-1)} \rvert \right). \$\$

The debiasing step refits the regression model on the active set
selected by the penalized estimator using weighted least squares,
reducing shrinkage bias and improving finite-sample interpretability of
the coefficient estimates.

This function depends on auxiliary routines (not shown here), including
`q_num_NA()` for Monte Carlo integration and `add_to_window()` for
sliding-window averaging.

## See also

[`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html),
[`simdata`](https://philchalmers.github.io/mirt/reference/simdata.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- farlr_debias(
  n = nrow(resp),
  resp = resp,
  parTab = parTab,
  K_hat = K_hat,
  p = p,
  lambda_all = seq(0.001, 0.1, length.out = 10),
  theta_est_irt.mean = theta_mean,
  theta_est_irt.se = theta_se,
  resp_rep = resp_rep,
  Z.em = Z.em,
  main = 1,
  verbose = TRUE
)
} # }
```
