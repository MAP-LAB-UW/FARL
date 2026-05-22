# FARLR EM-M Algorithm for Latent Regression with Regularization

Fits a FARLR-style latent regression model using an EM-M algorithm. In
the E-step, latent traits are sampled from a normal approximation based
on IRT estimates (`theta_est_irt.mean`, `theta_est_irt.se`). In the
M-step, regression coefficients are updated using weighted `glmnet`
(LASSO) with a debiasing refit via weighted least squares. A sliding
window averaging scheme is used to stabilize coefficient updates across
iterations. The tuning parameter `lambda` is selected by minimizing a
BIC-like criterion over `lambda_all`.

## Usage

``` r
farlr_emm(
  n,
  resp,
  parTab,
  K_hat,
  p,
  lambda_all,
  delta.criteria = 0.001,
  iter.max = 500,
  n_sam = 50,
  window.size = 50,
  theta_est_irt.mean,
  theta_est_irt.se,
  resp_rep,
  Z.em,
  Uupdate = NA,
  hatU = NA,
  Fan = NA,
  main,
  verbose = TRUE
)
```

## Arguments

- n:

  Integer. Sample size (number of persons).

- resp:

  Matrix. Observed item responses of dimension `n x J`.

- parTab:

  Data frame. Item parameter table containing at least `slope`,
  `difficulty`, and `guessin`.

- K_hat:

  Integer. Number of latent factors/components included in `Z.em`.

- p:

  Integer. Number of covariates/predictors included in `Z.em`.

- lambda_all:

  Numeric vector. Candidate regularization parameters passed to
  `glmnet`.

- delta.criteria:

  Numeric. Convergence tolerance for parameter updates. Default is
  `1e-3`.

- iter.max:

  Integer. Maximum number of EM iterations for each `lambda`. Default is
  `500`.

- n_sam:

  Integer. Number of Monte Carlo samples per subject used in the E-step.
  Default is `50`.

- window.size:

  Integer. Sliding window size used to average coefficient updates
  across iterations. Default is `50`.

- theta_est_irt.mean:

  Numeric vector of length `n`. IRT-based posterior mean estimates of
  latent trait \\\theta\\.

- theta_est_irt.se:

  Numeric vector of length `n`. IRT-based posterior standard error
  estimates of \\\theta\\.

- resp_rep:

  Matrix. Replicated/expanded responses used for Monte Carlo
  computations in the E-step (see `q_num_NA`).

- Z.em:

  Matrix. Design matrix used in the regression step, typically of
  dimension `n x (K_hat + p)`.

- main:

  Integer index (or indices). Predictor(s) to be treated as unpenalized
  via `penalty.factor` (set to 0). All other predictors are penalized
  unless already unpenalized in the first `K_hat` columns.

- verbose:

  Logical. If `TRUE`, prints progress messages and a progress bar.
  Default is `TRUE`.

## Value

A list with the following elements:

- `coefficients`:

  Estimated regression coefficients at the selected `lambda` (length
  `K_hat + p`).

- `sigma`:

  Estimated residual standard deviation.

- `LogLik`:

  BIC-like objective value corresponding to the selected `lambda` (named
  `LogLik` for compatibility).

- `minBIC`:

  Index of `lambda_all` achieving the minimum BIC criterion.

- `Convergence`:

  A character flag indicating convergence status.

## Details

For each `lambda` in `lambda_all`, the algorithm iterates until `delta`
falls below `delta.criteria` or `iter.max` is reached. The maximum
change in coefficient and residual scale estimates is monitored:
\\\delta = \max(\|\sigma^{(t)}-\sigma^{(t-1)}\|, \max_j
\|\beta_j^{(t)}-\beta_j^{(t-1)}\|)\\.

This function relies on helper functions (not shown here), including
`q_num_NA()` and `add_to_window()`.

## See also

[`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html),
[`simdata`](https://philchalmers.github.io/mirt/reference/simdata.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- farlr_emm(
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
